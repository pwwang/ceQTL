source("{{biopipen_dir}}/utils/misc.R")

library(snpStats)

exfile <- {{in.exfile | r}}
genofile <- {{in.genofile | r}}
snpgenefile <- {{out.snpgene | r}}
tftargetfile <- {{out.tftarget | r}}

log_info("Reading simulated expression data ...")
sim <- readRDS(exfile)
if (attr(sim, 'simulation_tool') == 'ESCO') {
    log_sim <- logNormCounts(sim, assay.type = "TrueCounts")
    exprs <- assays(log_sim)$logcounts
    rm(log_sim)
} else if (attr(sim, 'simulation_tool') == 'RUVcorr') {
    exprs <- t(sim$Y)
} else {
    stop("Unknown simulation tool.")
}
genes <- colnames(exprs)
rm(sim)
rm(exprs)

log_info("Reading genotype data ...")
# geno <- read.plink(file.path(gtdir, "sim_snps.bed"), file.path(gtdir, "sim_snps.bim"), file.path(gtdir, "sim_snps.fam"))
# geno <- as(geno$genotypes, "numeric")
geno <- read.table(
    genofile,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = 1)
snps <- colnames(geno)
rm(geno)

log_info("Extracting SNP-gene pairs ...")
snpgene <- data.frame(
    SNP = snps,
    ANNO = "",
    Gene = rep(paste(genes, collapse = "\t"), length(snps))
)
write.table(snpgene, snpgenefile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

log_info("Extracting TF-target pairs ...")
tftarget <- data.frame(
    TF = genes,
    ANNO = "",
    Target = rep(paste(genes, collapse = "\t"), length(genes))
)
write.table(tftarget, tftargetfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
