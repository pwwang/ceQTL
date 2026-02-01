#!/usr/bin/env fish

set ceqtl_genevar_pairs ./ceQTL-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set ceqtl_genevar_pairs_cov ./ceQTL-cov-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set goldfile ./eQTL-output/MatrixEQTL/cisqtls.txt
set goldfile_cov ./eQTL-cov-output/MatrixEQTL/cisqtls.txt
set gold_cutoff $argv[1]

poetry run -C ../../ceqtl -- \
    ceqtl roc \
        --name ceQTL_ROC_$gold_cutoff \
        --Input.in.infile $ceqtl_genevar_pairs \
        --Input.in.goldfile $goldfile \
        --Input.envs.in_ids SNP Target \
        --Input.envs.in_score MetaPval \
        --Input.envs.in_dir - \
        --Input.envs.gold_ids SNP gene \
        --Input.envs.gold_dir - \
        --Input.envs.gold_score FDR \
        --Input.envs.gold_cutoff $gold_cutoff

poetry run -C ../../ceqtl -- \
    ceqtl roc \
        --name ceQTL_ROC-cov_$gold_cutoff \
        --Input.in.infile $ceqtl_genevar_pairs_cov \
        --Input.in.goldfile $goldfile_cov \
        --Input.envs.in_ids SNP Target \
        --Input.envs.in_score MetaPval \
        --Input.envs.in_dir - \
        --Input.envs.gold_ids SNP gene \
        --Input.envs.gold_dir - \
        --Input.envs.gold_score FDR \
        --Input.envs.gold_cutoff $gold_cutoff
