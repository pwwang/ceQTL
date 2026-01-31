# ceQTL: A co-expression QTL model to detect a variant that affects transcription factor binding and its target regulation

Expression quantitative trait locus (eQTL) mapping is used to identify the functional link between a genomic variant, and a gene’s expression. A significant eQTL association does not mean causal relationship or mechanism and further investigation is needed to understand how a SNP impacts gene expression. One of the most plausible explanations for eQTL is that a genomic variant affects transcription factor (TF) binding and thus impacts its regulation on target genes (TGs). However, current eQTL does not provide information on the TF, and how its regulation is mediated by the SNP’s genotypes. Here, we propose a new method called differential co-expression QTL (ceQTL) among different alleles using Chow statistics to specifically detect eQTLs that are bound by a particular TF. We start with building a trio of TF, its TG, and related SNP and then test the significant coefficient difference among different genotypes of the SNP.  We applied this ceQTL model to simulated data and the lung tissue datasets from the GTEx project. The simulated data results showed that the model was robust to detect true ceQTLs at variable sample sizes and different minor allele frequencies as measured by Area Under the Curve (AUC). Our tool also performed a TF binding affinity analysis to add another layer of evidence for functional interpretation. In summary, ceQTL analysis provides a more interpretable and biological insight into the mechanism of eQTL and transcriptomic regulation, which would help us better understand how genomic variants affect phenotypes and diseases.

## Source Code

See the `ceqtl/` folder for the source code of the ceQTL model.

## Real Data

See the `real/` folder for the real data analysis using GTEx lung tissue data.
