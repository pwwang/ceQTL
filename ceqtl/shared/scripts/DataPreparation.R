source("{{biopipen_dir}}/utils/misc.R")

library(glue)
library(parallel)

genofile <- {{in.geno | r}}
exprfile <- {{in.expr | r}}
covfile <- {{in.cov | r}}
sgfile <- {{in.genesnp | r}}
tftfile <- {{in.tftarget | r}}
outdir <- {{out.outdir | r}}
ncores <- {{envs.ncores | r}}
nchunks <- {{envs.nchunks | r}}
transpose_geno <- {{envs.transpose_geno | r}}
transpose_expr <- {{envs.transpose_expr | r}}
transpose_cov <- {{envs.transpose_cov | r}}
set.seed(1234)

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
    close(con)
    return(gmt)
}

log_info("Reading genotype data ...")
geno <- read.table(
    genofile,
    row.names = 1,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
if (transpose_geno) { geno <- t(geno) }
allsnps <- colnames(geno)

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

samples <- intersect(rownames(geno), rownames(expr))
covdata <- NULL
cov <- NULL

if (is.null(covfile)) {
    log_info("- Working on {length(samples)} common samples between genotype and expression data.")
    geno <- geno[samples, , drop = FALSE]
    expr <- expr[samples, , drop = FALSE]
} else {
    log_info("Reading covariate data ...")
    covdata <- read.table(
        covfile,
        row.names = 1,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (transpose_cov) { covdata <- t(covdata) }
    ge_nsamples <- length(samples)
    samples <- intersect(samples, rownames(covdata))
    log_info("- Working on {length(samples)} common samples between genotype, expression and covariate data.")
    log_info("- Number of common samples between genotype and expression data: {ge_nsamples}")
    geno <- geno[samples, , drop = FALSE]
    expr <- expr[samples, , drop = FALSE]
    covdata <- covdata[samples, , drop = FALSE]
    cov <- paste0(bQuote(colnames(covdata)), collapse = " + ")
}

log_info("Reading SNP-gene pairs ...")
genesnp_pairs <- readGMT(sgfile)

log_info("Reading TF-target pairs ...")
tftarget_pairs <- readGMT(tftfile)

log_info("QC TF-target pairs ...")
nonexist_tfs <- c()
notarget_tfs <- c()
for (tf in names(tftarget_pairs)) {
    if (!tf %in% allgenes) {
        # log_warn("- TF {tf} does not exist in the expression data and will be ignored.")
        tftarget_pairs[[tf]] <- NULL
        nonexist_tfs <- c(nonexist_tfs, tf)
    } else {
        targets <- tftarget_pairs[[tf]]
        targets <- intersect(targets, allgenes)
        if (length(targets) == 0) {
            # log_warn("- TF {tf} has no targets in the expression data and will be ignored.")
            tftarget_pairs[[tf]] <- NULL
            notarget_tfs <- c(notarget_tfs, tf)
        }
        tftarget_pairs[[tf]] <- unique(targets)
    }
}
if (length(nonexist_tfs) > 0) {
    nonexist_tf_file <- file.path(outdir, "nonexist_tfs.txt")
    writeLines(nonexist_tfs, nonexist_tf_file)
    nonexist_tf_str <- if (length(nonexist_tfs) > 3) {
        paste0(paste(head(nonexist_tfs, 3), collapse = ", "), ", ...")
    } else {
        paste0(nonexist_tfs, collapse = ", ")
    }
    log_warn("- TFs do not exist in the expression data and will be ignored: ")
    log_warn("  {nonexist_tf_str} (n={length(nonexist_tfs)})")
    if (length(nonexist_tfs) > 3) {
        log_warn("  Full list {nonexist_tf_file}")
    }
}
if (length(notarget_tfs) > 0) {
    notarget_tf_file <- file.path(outdir, "notarget_tfs.txt")
    writeLines(notarget_tfs, notarget_tf_file)
    notarget_tf_str <- if (length(notarget_tfs) > 3) {
        paste0(paste(head(notarget_tfs, 3), collapse = ", "), ", ...")
    } else {
        paste0(notarget_tfs, collapse = ", ")
    }
    log_warn("- TFs have no targets in the expression data and will be ignored: ")
    log_warn("  {notarget_tf_str} (n={length(notarget_tfs)})")
    if (length(notarget_tfs) > 3) {
        log_warn("  Full list: {notarget_tf_file}")
    }
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

allgenes <- unique(c(names(tftarget_pairs), names(targettf_pairs)))
expr <- expr[, allgenes, drop = FALSE]

log_info("QC Gene-SNP pairs ...")
nonexist_snpgenes <- c()
nosnp_genes <- c()
for (gene in names(genesnp_pairs)) {
    if (!gene %in% allgenes) {
        # log_warn("- Gene {gene} does not exist in the expression data and will be ignored.")
        genesnp_pairs[[gene]] <- NULL
        nonexist_snpgenes <- c(nonexist_snpgenes, gene)
    } else {
        snps <- genesnp_pairs[[gene]]
        snps <- intersect(snps, allsnps)
        if (length(snps) == 0) {
            # log_warn("- Gene {gene} has no SNPs in the genotype data and will be ignored.")
            genesnp_pairs[[gene]] <- NULL
            nosnp_genes <- c(nosnp_genes, gene)
        }
        genesnp_pairs[[gene]] <- unique(snps)
    }
}
if (length(nonexist_snpgenes) > 0) {
    nonexist_snpgene_file <- file.path(outdir, "nonexist_snpgenes.txt")
    writeLines(nonexist_snpgenes, nonexist_snpgene_file)
    nonexist_snpgene_str <- if (length(nonexist_snpgenes) > 3) {
        paste0(paste(head(nonexist_snpgenes, 3), collapse = ", "), ", ...")
    } else {
        paste0(nonexist_snpgenes, collapse = ", ")
    }
    log_warn("- Genes do not exist in the expression data and will be ignored: ")
    log_warn("  {nonexist_snpgene_str} (n={length(nonexist_snpgenes)})")
    if (length(nonexist_snpgenes) > 3) {
        log_warn("  Full list: {nonexist_snpgene_file}")
    }
}
if (length(nosnp_genes) > 0) {
    nosnp_gene_file <- file.path(outdir, "nosnp_genes.txt")
    writeLines(nosnp_genes, nosnp_gene_file)
    nosnp_gene_str <- if (length(nosnp_genes) > 3) {
        paste0(paste(head(nosnp_genes, 3), collapse = ", "), ", ...")
    } else {
        paste0(nosnp_genes, collapse = ", ")
    }
    log_warn("- Genes have no SNPs in the genotype data and will be ignored: ")
    log_warn("  {nosnp_gene_str} (n={length(nosnp_genes)})")
    if (length(nosnp_genes) > 3) {
        log_warn("  Full list: {nosnp_gene_file}")
    }
}

# Reverse the gene-SNP pair mapping
snpgene_pairs <- list()
for (gene in names(genesnp_pairs)) {
    for (snp in genesnp_pairs[[gene]]) {
        if (is.null(snpgene_pairs[[snp]])) {
            snpgene_pairs[[snp]] <- gene
        } else {
            snpgene_pairs[[snp]] <- unique(c(snpgene_pairs[[snp]], gene))
        }
    }
}

log_info("Splitting genes into chunks ...")
clusters <- sample(1:nchunks, length(genesnp_pairs), replace = TRUE)
km <- list(cluster = clusters)
names(km$cluster) <- names(genesnp_pairs)
clustered = list(km = km, clusters = sort(unique(clusters)))

process_one_cluster <- function(cluster) {
    log_info("- Processing cluster {cluster}/{length(clustered$clusters)} ...")
    chunk_genges <- names(genesnp_pairs)[clustered$km$cluster == cluster]
    chunk_snps <- unique(unlist(lapply(chunk_genges, function(gene) genesnp_pairs[[gene]])))
    chunk_geno <- geno[, chunk_snps , drop = FALSE]
    chunk_dir <- file.path(
        outdir,
        paste0(tools::file_path_sans_ext(basename(genofile)), ".chunk-", cluster)
    )
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
    chunk_allgenes <- c()
    chow_lines <- c(
        # Target ~ TF + Covariates
        paste("Group", "Formula", "SNP", "TF", "Target", sep = "\t"),
        sapply(chunk_snps, function(snp) {
            targets <- snpgene_pairs[[snp]]
            chunk_allgenes <<- union(chunk_allgenes, targets)
            sapply(targets, function(target) {
                tfs <- unique(targettf_pairs[[target]])
                chunk_allgenes <<- union(chunk_allgenes, tfs)
                sapply(tfs, function(tf) {
                    chow_formula <- paste0(bQuote(target), " ~ ", bQuote(tf))
                    if (!is.null(cov)) {
                        chow_formula <- paste0(chow_formula, " + ", cov)
                    }
                    paste0(snp, "\t", chow_formula, "\t", snp, "\t", tf, "\t", target)
                })
            })
        })
    )
    interactionlm_lines <- c(
        # Target ~ TF*SNP + Covariates
        paste("SNP", "TF", "Target", "Formula", sep = "\t"),
        sapply(chunk_snps, function(snp) {
            targets <- snpgene_pairs[[snp]]
            # chunk_allgenes <<- union(chunk_allgenes, targets)
            sapply(targets, function(target) {
                tfs <- unique(targettf_pairs[[target]])
                # chunk_allgenes <<- union(chunk_allgenes, tfs)
                sapply(tfs, function(tf) {
                    interactionlm_formula <- paste0(
                        bQuote(target), " ~ ", bQuote(tf), " * ", bQuote(snp)
                    )
                    if (!is.null(cov)) {
                        interactionlm_formula <- paste0(interactionlm_formula, " + ", cov)
                    }
                    paste0(snp, "\t", tf, "\t", target, "\t", interactionlm_formula)
                })
            })
        })
    )
    la_lines <- sapply(chunk_snps, function(snp) {
        targets <- snpgene_pairs[[snp]]
        # chunk_allgenes <<- union(chunk_allgenes, targets)
        sapply(targets, function(target) {
            tfs <- unique(targettf_pairs[[target]])
            # chunk_allgenes <<- union(chunk_allgenes, tfs)
            sapply(tfs, function(tf) paste0(snp, "\t", tf, "\t", target))
        })
    })

    # Save the chow formula
    writeLines(unlist(chow_lines), file.path(chunk_dir, "formula-chow.txt"))
    # Save the interactionlm formula
    writeLines(unlist(interactionlm_lines), file.path(chunk_dir, "formula-interactionlm.txt"))
    # Save the liquid association formula
    writeLines(unlist(la_lines), file.path(chunk_dir, "formula-la.txt"))

    # subset the expression data
    chunk_expr <- expr[, unique(chunk_allgenes), drop = FALSE]
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

tryCatch({
    x <- mclapply(clustered$clusters, process_one_cluster, mc.cores = ncores)
}, warning = function(w) {
    if (grepl("errors", w)) {
        log_error(x[[1]])
        stop(w)
    } else {
        log_warn(w)
    }
})
