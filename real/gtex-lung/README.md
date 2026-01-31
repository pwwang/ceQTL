# ceQTLs on GTex Lung tissue

This folder contains the ceQTLs identified in the GTex Lung tissue and downstream analysis.

## Steps

### Data Preparation

1. Quality control of the genotype data, keep the SNPs in the promoter regions (up/down 2500bp). Also generate the coordinate file for the SNPs and Gene-SNP pairs.
2. Filter the gene expression data, keep the genes in the TF-Gene network.

### Calling ceQTLs

3. Calling w/o covariates

### Calling eQTLs

Only cis-eQTLs are called in this step (up/down 2500bp).

4. Calling w/o covariates

### Manhatton Plot

5. Plotting the Manhatton plot for the ceQTLs w/o covariates

### ROC against eQTLs

6. ROC for the ceQTLs w/o covariates

### QQ Plots

7. QQ plot for the ceQTLs w/o covariates

## Results

The results in directory with "-cov-output" are the counterparts of the results in the directory without "-cov-output".

The pvalues or adjusted pvalues (FDR method used) for ceQTL trios are the results directly from the model.
The (meta) pvalues or adjusted (meta) pvalues for ceQTL variants are combined using fisher's method.

- ceQTLs:
  - Trios (SNP-TF-Target): ceQTL-output/CombineAndAdjustPValues/ceqtl-trios.txt
  - SNP-Gene pairs: ceQTL-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
  - SNPs: ceQTL-output/CombinePValuesForVariants/ceqtl-variants.txt
- eQTLs:
  - eQTL-output/MatrixEQTL/cisqtls.txt
- Manhatton plots:
  - Use ceQTL trio pvalues: ManhattanPlots-output/ManhattanTrios/ceqtl-trios.manhattan.png
  - Use ceQTL trio adjusted pvalues: ManhattanPlots_Padj-output/ManhattanTrios/ceqtl-trios.manhattan.png
- ROC, use eQTL as gold standard (positives: FDR < 1e-10):
  - ceQTL_ROC-output/ROC/ceqtl-vargenes.4roc.roc.png
- QQ plots:
  - For trios + Pval: QQPlot_trios_Pval-output/QQPlot/ceqtl-trios.qq.png
  - For trios + Padj: QQPlot_trios_Padj-output/QQPlot/ceqtl-trios.qq.png
  - For variants + Pval: QQPlot_vars_Pval-output/QQPlot/ceqtl-variants.qq.png
  - For variants + Padj: QQPlot_vars_Padj-output/QQPlot/ceqtl-variants.qq.png
