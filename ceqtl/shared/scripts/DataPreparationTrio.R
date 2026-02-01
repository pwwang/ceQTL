library(biopipen.utils)
library(parallel)

genofile <- {{in.geno | r}}
exprfile <- {{in.expr | r}}
covfile <- {{in.cov | r}}
triofile <- {{in.triofile | r}}
outdir <- {{out.outdir | r}}
transpose_geno <- {{envs.transpose_geno | r}}
transpose_expr <- {{envs.transpose_expr | r}}
transpose_cov <- {{envs.transpose_cov | r}}
set.seed(1234)

log <- get_logger()

log$info("Reading genotype data ...")
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

log$info("Reading expression data ...")
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
    log$info("- Working on {length(samples)} common samples between genotype and expression data.")
    geno <- geno[samples, , drop = FALSE]
    expr <- expr[samples, , drop = FALSE]
} else {
    log$info("Reading covariate data ...")
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
    log$info("- Working on {length(samples)} common samples between genotype, expression and covariate data.")
    log$info("- Number of common samples between genotype and expression data: {ge_nsamples}")
    geno <- geno[samples, , drop = FALSE]
    expr <- expr[samples, , drop = FALSE]
    covdata <- covdata[samples, , drop = FALSE]
    cov <- paste0(sapply(colnames(covdata), bQuote), collapse = " + ")
}

log$info("Reading trio data ...")
trios <- read.table(
    triofile,
    row.names = NULL,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
cnames <- tolower(colnames(trios))
tfcol <- match("tf", cnames)
snpcol <- match("snp", cnames)
targetcol <- setdiff(1:3, c(tfcol, snpcol))
trios <- trios[, c(tfcol, snpcol, targetcol)]
colnames(trios) <- c("tf", "snp", "target")

chunk_dir <- file.path(
    outdir,
    paste0(tools::file_path_sans_ext(basename(genofile)), ".chunk-", 1)
)
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

chunk_geno <- geno[, intersect(allsnps, unique(trios$snp)), drop = FALSE]
# Save genotype data
write.table(
    chunk_geno,
    file.path(chunk_dir, "genotype.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
)

chunk_expr <- expr[, intersect(allgenes, unique(c(trios$tf, trios$target))), drop = FALSE]
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

chow_df <- data.frame(
    Group = character(),
    Formula = character(),
    SNP = character(),
    TF = character(),
    Target = character()
)
intlm_df <- data.frame(
    SNP = character(),
    TF = character(),
    Target = character(),
    Formula = character()
)
la_df <- data.frame(
    SNP = character(),
    TF = character(),
    Target = character()
)

chunk_genes <- colnames(chunk_expr)
chunk_snps <- colnames(chunk_geno)

apply(trios, 1, function(row) {
    tf <- row[1]
    snp <- row[2]
    target <- row[3]
    if (tf %in% chunk_genes && snp %in% chunk_snps && target %in% chunk_genes) {
        chow_df <<- rbind(chow_df, data.frame(
            Group = snp,
            Formula = if (is.null(cov)) {
                paste0(bQuote(target), " ~ ", bQuote(tf))
            } else {
                paste0(bQuote(target), " ~ ", bQuote(tf), " + ", cov)
            },
            SNP = snp,
            TF = tf,
            Target = target
        ))
        intlm_df <<- rbind(intlm_df, data.frame(
            SNP = snp,
            TF = tf,
            Target = target,
            Formula = if (is.null(cov)) {
                paste0(bQuote(target), " ~ ", bQuote(tf), " * ", bQuote(snp))
            } else {
                paste0(bQuote(target), " ~ ", bQuote(tf), " * ", bQuote(snp), " + ", cov)
            }
        ))
        la_df <<- rbind(la_df, data.frame(
            SNP = snp,
            TF = tf,
            Target = target
        ))
    }
})

log$info("Writing formula files ...")
write.table(chow_df, file.path(chunk_dir, "formula-chow.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(intlm_df, file.path(chunk_dir, "formula-interactionlm.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(la_df, file.path(chunk_dir, "formula-la.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
