source("{{biopipen_dir}}/utils/misc.R")

genofile <- {{in.geno | r}}
exprfile <- {{in.expr | r}}
covfile <- {{in.cov | r}}
sgfile <- {{in.snpgene | r}}
tftfile <- {{in.tftarget | r}}
outdir <- {{out.outdir | r}}
nchunks <- {{envs.nchunks | r}}
transpose_geno <- {{envs.transpose_geno | r}}
transpose_expr <- {{envs.transpose_expr | r}}

readGMT <- function(gmtfile) {
    gmt <- list()
    con <- file(gmtfile, "r")
    while (TRUE) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) {
            break
        }
        items <- strsplit(trimws(line), "\t")[[1]]
        gmt[[items[1]]] <- items[3:length(items)]
    }
    return(gmt)
}

log_info("Reading genotype data ...")
geno <- read.table(
    genofile,
    row.names = 1
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
if (transpose_geno) { geno <- t(geno) }

log_info("Reading expression data ...")
expr <- read.table(
    exprfile,
    row.names = 1,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
if (transpose_expr) { expr <- t(expr) }
allgenes <- colnames(expr)

if (nrow(geno) != nrow(expr)) {
    stop("Number of samples in genotype and expression data do not match.")
}

if (!isTRUE(all.equal(rownames(geno), rownames(expr)))) {
    stop("Sample IDs in genotype and expression data do not match.")
}

covdata <- NULL
cov <- NULL
if (!is.null(covfile)) {
    log_info("Reading covariate data ...")
    covdata <- read.table(
        covfile,
        row.names = 1,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (nrow(geno) != nrow(covdata)) {
        stop("Number of samples in genotype and covariate data do not match.")
    }
    if (!isTRUE(all.equal(rownames(geno), rownames(covdata)))) {
        stop("Sample IDs in genotype and covariate data do not match.")
    }
    cov <- paste0(bQuote(colnames(covdata)), collapse = " + ")
}

log_info("Reading SNP-gene pairs ...")
snpgene_pairs <- readGMT(sgfile)

log_info("Reading TF-target pairs ...")
tftarget_pairs <- readGMT(tftfile)
nonexisting_tfs <- setdiff(names(tftarget_pairs), allgenes)
if (length(nonexisting_tfs) > 0) {
    log_warn("- The following TFs do not exist in the expression data and will be ignored:")
    log_warn("- {paste(head(nonexisting_tfs), collapse = ', ')} (total={length(nonexisting_tfs)}) ...")
    tftarget_pairs <- tftarget_pairs[setdiff(names(tftarget_pairs), nonexisting_tfs)]
}
# Reverse the TF-target pair mapping
targettf_pairs <- list()
for (tf in names(tftarget_pairs)) {
    for (target in tftarget_pairs[[tf]]) {
        if (is.null(targettf_pairs[[target]])) {
            targettf_pairs[[target]] <- tf
        } else {
            targettf_pairs[[target]] <- unique(c(targettf_pairs[[target]], tf))
        }
    }
}
nonexisting_targets <- setdiff(names(targettf_pairs), allgenes)
if (length(nonexisting_targets) > 0) {
    log_warn("- The following targets do not exist in the expression data and will be ignored:")
    log_warn("- {paste(head(nonexisting_targets), collapse = ', ')} (total={length(nonexisting_targets)}) ...")
    targettf_pairs <- targettf_pairs[setdiff(names(targettf_pairs), nonexisting_targets)]
}
allgenes <- c(names(tftarget_pairs), names(targettf_pairs))
expr <- expr[, allgenes, drop = FALSE]

log_info("Clustering SNPs and splitting into chunks ...")
# Used for clustering
snpgene_indicators <- list()
for (snp in names(snpgene_pairs)) {
    # remove genes without TFs
    snpgene_pairs[[snp]] <- unique(intersect(snpgene_pairs[[snp]], names(targettf_pairs)))
    tfs <- unlist(unname(targettf_pairs[snpgene_pairs[[snp]]]))
    snpgene_indicators[[snp]] <- allgenes %in% c(snpgene_pairs[[snp]], tfs)
}
snpgene_indicators <- do_call(rbind, snpgene_indicators)
km <- kmeans(snpgene_indicators, nchunks)
clusters <- unique(km$cluster)

for (cluster in clusters) {
    log_info("- Processing cluster {cluster}/{length(clusters)} ...")
    chunk_snps <- names(km$cluster[km$cluster == cluster])
    chunk_geno <- geno[chunk_snps, , drop = FALSE]
    chunk_dir <- file.path(outdir, paste0("chunk-", cluster))
    dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

    # Save genotype data
    write.table(
        chunk_geno,
        file.path(chunk_dir, "genotype.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
    )

    # Compose the formula for each SNP
    # Target ~ TF + Covariates
    chow_lines <- paste("Group", "Formula", "SNP", "TF", "Target", sep = "\t")
    # Target ~ TF*SNP + Covariates
    interactionlm_lines <- paste("SNP", "TF", "Target", "Formula", sep = "\t")
    chunk_allgenes <- c()
    for (snp in chunk_snps) {
        targets <- snpgene_pairs[[snp]]
        chunk_allgenes <- c(chunk_allgenes, targets)
        for (target in targets) {
            tfs <- unique(targettf_pairs[[target]])
            chunk_allgenes <- c(chunk_allgenes, tfs)
            for (tf in tfs) {
                chow_formula <- paste0(bQuote(target), " ~ ", bQutoe(tf))
                interactionlm_formula <- paste0(
                    bQuote(target), " ~ ", bQuote(tf), " * ", bQuote(snp)
                )
                if (!is.null(cov)) {
                    chow_formula <- paste0(chow_formula, " + ", cov)
                    interactionlm_formula <- paste0(interactionlm_formula, " + ", cov)
                }
                chow_lines <- c(
                    chow_lines,
                    paste0(snp, "\t", chow_formula, "\t", snp, "\t", tf, "\t", target)
                )
                interactionlm_lines <- c(
                    interactionlm_lines,
                    paste0(snp, "\t", tf, "\t", target, "\t", interactionlm_formula)
                )
            }
        }
    }
    # Save the chow formula
    writeLines(chow_lines, file.path(chunk_dir, "formula-chow.txt"))
    # Save the interactionlm formula
    writeLines(interactionlm_lines, file.path(chunk_dir, "formula-interactionlm.txt"))

    # subset the expression data
    chunk_expr <- expr[unique(chunk_allgenes), , drop = FALSE]
    # Save the expression data
    write.table(
        chunk_expr,
        file.path(chunk_dir, "expression-nocov.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
    )

    # add covariates
    if (!is.null(covdata)) {
        chunk_expr <- cbind(chunk_expr, covdata)
    }

    # Save the expression data
    write.table(
        chunk_expr,
        file.path(chunk_dir, "expression.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
    )
}
