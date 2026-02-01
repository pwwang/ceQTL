#!/usr/bin/env fish

set ceqtl_trios ./ceQTL-output/CombineAndAdjustPValues/ceqtl-trios.txt
set ceqtl_vars ./ceQTL-output/CombinePValuesForVariants/ceqtl-variants.txt
set ceqtl_trios_cov ./ceQTL-cov-output/CombineAndAdjustPValues/ceqtl-trios.txt
set ceqtl_vars_cov ./ceQTL-cov-output/CombinePValuesForVariants/ceqtl-variants.txt
set bedfile ./ProcessGT-output/VariantCoords/variant_coords.bed
set chroms 'chr1-22'
set signif '1e-10, 1e-8'

# without cov
poetry run -C ../../ceqtl -- \
	ceqtl manh --name ManhattanPlots_Pval \
	--ceqtl-trios $ceqtl_trios --ceqtl-vars $ceqtl_vars --bedfile $bedfile \
	--chroms $chroms \
	--signif $signif \
	--trio-pval-col Pval \
	--var-pval-col MetaPval \
	--signif 1e-10,1e-6

poetry run -C ../../ceqtl -- \
	ceqtl manh --name ManhattanPlots_Padj \
	--ceqtl-trios $ceqtl_trios --ceqtl-vars $ceqtl_vars --bedfile $bedfile \
	--chroms $chroms \
	--signif $signif \
	--trio-pval-col Padj \
	--var-pval-col MetaPadj \
	--signif 1e-10,1e-6

# with cov
poetry run -C ../../ceqtl -- \
	ceqtl manh --name ManhattanPlots_Pval-cov \
	--ceqtl-trios $ceqtl_trios_cov --ceqtl-vars $ceqtl_vars_cov --bedfile $bedfile \
	--chroms $chroms \
	--signif $signif \
	--trio-pval-col Pval \
	--var-pval-col MetaPval \
	--signif 1e-10,1e-6

poetry run -C ../../ceqtl -- \
	ceqtl manh --name ManhattanPlots_Padj-cov \
	--ceqtl-trios $ceqtl_trios_cov --ceqtl-vars $ceqtl_vars_cov --bedfile $bedfile \
	--chroms $chroms \
	--signif $signif \
	--trio-pval-col Padj \
	--var-pval-col MetaPadj \
	--signif 1e-10,1e-6
