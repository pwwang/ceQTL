#!/usr/bin/env fish

set ceqtl_genevar_pairs ./ceQTL-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set ceqtl_genevar_pairs_cov ./ceQTL-cov-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set goldfile ./eQTL-output/MatrixEQTL/cisqtls.txt
set goldfile_cov ./eQTL-cov-output/MatrixEQTL/cisqtls.txt
set trans " -log10"

poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot --name QQPlot_Against_eQTLs \
        --in.infile $ceqtl_genevar_pairs \
        --in.theorfile $goldfile \
        --envs.val_col MetaPval \
        --envs.theor_col p-value \
        --envs.trans $trans \
        --envs.theor_trans $trans \
        --envs.args.distribution custom \
        --envs.band.disabled


poetry run -C ../../ceqtl -- \
    pipen run plot QQPlot --name QQPlot_Against_eQTLs-cov \
        --in.infile $ceqtl_genevar_pairs_cov \
        --in.theorfile $goldfile_cov \
        --envs.val_col MetaPval \
        --envs.theor_col p-value \
        --envs.trans $trans \
        --envs.theor_trans $trans \
        --envs.args.distribution custom \
        --envs.band.disabled \
        --envs.line.distribution norm \
        --envs.line.disabled \
        --envs.ggs 'geom_abline(slope = 1, intercept = 0)'
