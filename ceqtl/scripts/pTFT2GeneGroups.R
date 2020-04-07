{{"__init__.R" | rimport}}
library(reshape2)
infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}

# infile is like:
# TF6	Gene1
# TF5	Gene120
# TF7	Gene76
# TF18	Gene155
# TF11	Gene931
# TF12	Gene478
# TF3	Gene142

# turn it into:
#       Gene1 Gene2
# TF1   TF  NA
# TF2   NA  TF
# Gene1 Gene  NA
# Gene2 NA  Gene
# ...

indata = read.table.inopts(infile, list(cnames=FALSE, rnames=FALSE))
colnames(indata) = c("TF", "Target")

wide = dcast(indata, TF ~ Target, value.var = "Target")
tfs = wide[, 1, drop=TRUE]
rownames(wide) = tfs
wide = wide[, -1]
wide[!is.na(wide)] = "TF"
genes = colnames(wide)

# attach the diagnal to wide
ncols = ncol(wide) # number of genes
dgmat = matrix(NA, ncol = ncols, nrow = ncols)
diag(dgmat) = "Target"
rownames(dgmat) = colnames(dgmat) = genes

# merge same genes (rows)
for (tf in tfs) {
    if (!tf %in% genes) {
        next
    }
    wide[tf, is.na(wide[tf, ])] = dgmat[, is.na(wide[tf, ])]
}
wide = rbind(wide, dgmat[setdiff(genes, tfs), , drop= FALSE])

write.table(wide, outfile, sep="\t", quote=FALSE)
