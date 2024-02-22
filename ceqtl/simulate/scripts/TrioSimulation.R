source("{{biopipen_dir}}/utils/misc.R")
library(dcanr)
library(scuttle)
library(doRNG)
library(doParallel)
library(snpStats)
library(rlang)
library(dplyr)

exfile <- {{in.exfile | r}}
gtdir <- {{in.gtdir | r}}
exprfile <- {{out.exprfile | r}}
genofile <- {{out.genofile | r}}
triofile <- {{out.triofile | r}}
seed <- {{envs.seed | r}}
padj <- {{envs.padj | r}}
keep <- {{envs.keep | r}}
beta <- {{envs.beta | r}}
ncores <- {{envs.ncores | r}}
cor_method <- {{envs.cor_method | r}}
perm_batch <- {{envs.perm_batch | r}}

log_info("Setting seed and parallel backend ...")
set.seed(seed)
registerDoParallel(cores = ncores)
registerDoRNG(seed)

log_info("Reading simulated expression data ...")
#       Cell1 Cell2 Cell3 Cell4 Cell5
# Gene1     1     3     0     0     0
# Gene2     0     5     5     3     3
# Gene3     0     0     0     0     0
# Gene4     0     0     0     0     0
# Gene5     3     5     0     3     0
sim <- readRDS(exfile)
if (attr(sim, 'simulation_tool') == 'ESCO') {
    log_sim <- logNormCounts(sim, assay.type = "TrueCounts")
    exprs <- assays(log_sim)$logcounts
} else if (attr(sim, 'simulation_tool') == 'RUVcorr') {
    exprs <- t(sim$Y)
} else {
    stop("Unknown simulation tool.")
}
samples <- colnames(exprs)

log_info("Saving expression data ...")
write.table(t(exprs), exprfile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

log_info("Reading genotype data ...")
geno <- read.plink(file.path(gtdir, "sim_snps.bed"), file.path(gtdir, "sim_snps.bim"), file.path(gtdir, "sim_snps.fam"))
geno <- as(geno$genotypes, "numeric")
rownames(geno) <- samples
snps <- colnames(geno)
nsnps <- length(snps)

log_info("Saving genotype data ...")
write.table(geno, genofile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# extract trios for each snp
log_info("Extracting trios ...")

diffcoex.score <- function(snp, diffcoex.beta = beta, ...) {

    samples0 <- geno[, snp, drop = FALSE] == 0
    samples0 <- names(samples0[samples0,])
    expr0 <- exprs[, samples0, drop = FALSE]
    samples1 <- geno[, snp, drop = FALSE] == 1
    samples1 <- names(samples1[samples1,])
    expr1 <- exprs[, samples1, drop = FALSE]
    samples2 <- geno[, snp, drop = FALSE] == 2
    samples2 <- names(samples2[samples2,])
    expr2 <- exprs[, samples2, drop = FALSE]
    if (min(c(length(samples0), length(samples1), length(samples2))) < 3) {
        log_warn("  Less than 3 samples in one of the genotypes. Skipping ...")
        return(NULL)
    }
    rm(samples0, samples1, samples2)

    # calculate correlation
    r0 <- dcanr:::cor.pairs(t(expr0), cor.method = cor_method)
    r1 <- dcanr:::cor.pairs(t(expr1), cor.method = cor_method)
    r2 <- dcanr:::cor.pairs(t(expr2), cor.method = cor_method)

    C0 <- (sign(r0) * r0^2 + sign(r1) * r1^2 + sign(r2) * r2^2) / 3
    D <- abs(sign(r0) * r0^2 - C0) + abs(sign(r1) * r1^2 - C0) + abs(sign(r2) * r2^2 - C0)
    D <- sqrt(0.25 * D)
    D <- D^diffcoex.beta
    T.ovlap <- D %*% D + ncol(D) * D  #calc topological ovlap

    mins = matrix(rep(rowSums(D), ncol(D)), nrow = ncol(D))
    mins = pmin(mins, matrix(rep(colSums(D), each = ncol(D)), nrow = ncol(D)))
    T.ovlap = 1 - (T.ovlap/(mins + 1 - D))

    diag(T.ovlap) = 1

    #add run parameters as attributes
    attributes(T.ovlap) = c(
        attributes(T.ovlap),
        'cor.method' = cor_method,
        'diffcoex.beta' = diffcoex.beta,
        'call' = match.call()
    )

    return(1 - T.ovlap)
}

perm.test <- function(dcscores, snp, B = perm_batch) {
    obs = dcanr:::mat2vec(dcscores)

    #package requirements
    pckgs = c('dcanr')

    #perform permutation
    pvals = foreach(
        b = seq_len(B),
        .combine = function(...) {mapply(sum, ...)},
        .multicombine = TRUE,
        .inorder = FALSE,
        .packages = pckgs
    ) %dorng% {
        #shuffle condition and recalculate scores
        env = new.env()
        assign('snp', snp, envir = env)
        permsc = eval(attr(dcscores, 'call'), envir = env)
        permsc = dcanr:::mat2vec(permsc)

        #count elements greater than obs
        permsc = abs(permsc)
        permsc = permsc[!(is.na(permsc) || is.infinite(permsc))]
        permcounts = vapply(abs(obs), function(x) sum(permsc > x), 0)
        return(c(permcounts, length(permsc)))
    }

    #p-values
    N <- pvals[length(pvals)]
    pvals <- pvals[-(length(pvals))] / N
    # attributes(pvals) = attributes(obs)
    # pvals = dcanr:::vec2mat(pvals)
    # attr(pvals, 'dc.test') = 'permutation'
    # return(pvals)
    # Format into SNP,Gene1,Gene2,Pval
    gene_pairs <- as.data.frame(t(combn(attr(obs, 'feature.names'), 2)))
    colnames(gene_pairs) <- c('Gene1', 'Gene2')
    gene_pairs$SNP <- snp
    gene_pairs$Pval <- pvals
    gene_pairs[, c('SNP', 'Gene1', 'Gene2', 'Pval'), drop = FALSE]
}

do_one_snp <- function(i) {
    snp <- snps[i]
    log_info("- Processing SNP {i}/{nsnps}: {snp} ...")
    log_info("  Calculating differential co-expression scores ...")
    dcscores <- diffcoex.score(snp)

    if (!is.null(dcscores)) {
        log_info("  Calculating p-values ...")
        perm.test(dcscores, snp)
    }
}

trios <- do_call(rbind, lapply(seq_len(nsnps), do_one_snp))
trios$Padj <- p.adjust(trios$Pval, method = padj)
trios <- trios %>% filter(!!parse_expr(keep))

log_info("Saving trios ...")
write.table(trios, triofile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
