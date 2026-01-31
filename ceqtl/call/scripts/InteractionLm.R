source("{{biopipen_dir}}/utils/misc.R")

infile <- {{in.infile | r}}
groupfile <- {{in.groupfile | r}}
fmlfile <- {{in.fmlfile | r}}
outfile <- {{out.outfile | r}}
padj <- {{envs.padj | r}}

log_info("Reading input data ...")
indata <- read.table(
    infile,
    row.names = 1,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
groupdata <- read.table(
    groupfile,
    row.names = 1,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
fmldata <- read.table(
    fmlfile,
    row.names = NULL,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

indata <- cbind(indata, groupdata)

log_info("Running tests ...")
ncases <- nrow(fmldata)
results <- do_call(rbind, lapply(
    seq_len(ncases),
    function(i) {
		fmlrow <- fmldata[i, , drop=TRUE]
        if (i %% 100 == 0) {
            log_info("- {i} / {ncases} ...")
        }
        log_debug("  Running test for formula: {fmlrow$Formula}")

        lmfit <- lm(as.formula(fmlrow$Formula), data = indata)

        #ano <- anova(lmfit)
        #if (nrow(ano) < 4) {
        #    # `SPI1` ~ `SPI1` * `chr11_47376163_C_G_b38
        #    return(NULL)
        #}
        #tf_pval = ano[1, 5]
        #snp_pval = ano[2, 5]
        #intx_pval = ano[3, 5]
        ## Extract the coefficient
        #coefs <- coef(summary(lmfit))
        #tf_coef <- coefs[2, "Estimate"]
        #snp_coef <- coefs[3, "Estimate"]
        #intx_coef <- coefs[4, "Estimate"]

        ret <- summary(lmfit)$coefficients
        tf_pval = ret[2, 4]
        snp_pval = ret[3, 4]
        intx_pval = ret[nrow(ret), 4]
        tf_coef <- ret[2, 1]
        snp_coef <- ret[3, 1]
        intx_coef <- ret[nrow(ret), 1]

		fmlrow$Pval <- intx_pval
        fmlrow$Tf_Pval <- tf_pval
        fmlrow$Snp_Pval <- snp_pval
        fmlrow$Tf_Coef <- tf_coef
        fmlrow$Snp_Coef <- snp_coef
        fmlrow$Intx_Coef <- intx_coef
		fmlrow
    }
))
results <- as.data.frame(results)

if (padj != "none") {
    log_info("Adjusting p-values ...")
    results$Padj <- p.adjust(results$Pval, method = padj)
    results$Tf_Padj <- p.adjust(results$Tf_Pval, method = padj)
    results$Snp_Padj <- p.adjust(results$Snp_Pval, method = padj)
}

log_info("Writing output ...")
results <- apply(results, 2, as.character)
write.table(
    results,
    file = outfile,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)
