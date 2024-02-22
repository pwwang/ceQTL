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
        pval = anova(lmfit)[3, 5]
		fmlrow$Pval <- pval
		fmlrow
    }
))
results <- as.data.frame(results)

if (padj != "none") {
    log_info("Adjusting p-values ...")
    results$Padj <- p.adjust(results$Pval, method = padj)
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
