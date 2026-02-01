library(dplyr)
library(tidyr)
library(plotthis)

gtfile <- "ProcessGT-output/Plink2GTMat/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.txt"
gtmat <- read.table(gtfile, header=TRUE, row.names=1)

gtmat <- gtmat %>%
    rowwise() %>%
    mutate(n_refhom = sum(c_across(everything()) == 0),
           n_het = sum(c_across(everything()) == 1),
           n_althom = sum(c_across(everything()) == 2))

histdat <- data.frame(
    n_refhom = gtmat$n_refhom,
    n_het = gtmat$n_het,
    n_althom = gtmat$n_althom
)

histdat <- histdat %>%
    pivot_longer(cols = everything(),
                 names_to = "genotype",
                 values_to = "count")

print(histdata %>% group_by(genotype) %>%
          summarise(mean_count = mean(count),
                    min_count = min(count),
                    max_count = max(count),
                    median_count = median(count),
                    sd_count = sd(count)))

p <- Histogram(
    data = histdat,
    x = "count",
    binwidth = 10,
    group_by = "genotype")

outdir <- "Histogram-NSamples-output"
dir.create(outdir, showWarnings = FALSE)
png(file.path(outdir, "n_samples_histogram.png"), width=1800, height=1200, res=300)
print(p)
dev.off()
