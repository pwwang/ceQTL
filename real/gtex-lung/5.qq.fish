#!/usr/bin/env fish

set ceqtl_trios ./ceQTL-output/CombineAndAdjustPValues/ceqtl-trios.txt
set ceqtl_vars ./ceQTL-output/CombinePValuesForVariants/ceqtl-variants.txt
set ceqtl_trios_cov ./ceQTL-cov-output/CombineAndAdjustPValues/ceqtl-trios.txt
set ceqtl_vars_cov ./ceQTL-cov-output/CombinePValuesForVariants/ceqtl-variants.txt
set distribution '{"distribution": "exp"}'

poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_trios_Pval \
    --in.infile $ceqtl_trios \
    --envs.val_col Pval \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL trios: -log10(Pval)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_vars_Pval \
    --in.infile $ceqtl_vars \
    --envs.val_col MetaPval \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL variants: -log10(Combined Pval)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_trios_Padj \
    --in.infile $ceqtl_trios \
    --envs.val_col Padj \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL trios: -log10(Padj)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_vars_Padj \
    --in.infile $ceqtl_vars \
    --envs.val_col MetaPadj \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL variants: -log10(Combined Pvals)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_trios_Pval-cov \
    --in.infile $ceqtl_trios_cov \
    --envs.val_col Pval \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL trios: -log10(Pval)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_vars_Pval-cov \
    --in.infile $ceqtl_vars_cov \
    --envs.val_col MetaPval \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL variants: -log10(Combined Pval)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_trios_Padj-cov \
    --in.infile $ceqtl_trios_cov \
    --envs.val_col Padj \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL trios: -log10(Padj)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot \
    --name QQPlot_vars_Padj-cov \
    --in.infile $ceqtl_vars_cov \
    --envs.val_col MetaPadj \
    --envs.trans " -log10" \
    --envs.title "QQ Plot on ceQTL variants: -log10(Combined Pvals)" \
    --envs.band $distribution \
    --envs.line $distribution \
    --envs.point $distribution
