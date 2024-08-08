#!/usr/bin/env fish

set exprfile ProcessExpr-output/FilterGenes/Lung.v8.normalized_expression.matrix.txt
set genofile ProcessGT-output/Plink2GTMat/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.txt
set covfile raw/Lung.v8.covariates.txt
set ceqtl_trios ./ceQTL-output/CombineAndAdjustPValues/ceqtl-trios.txt
set ceqtl_trios_cov ./ceQTL-cov-output/CombineAndAdjustPValues/ceqtl-trios.txt
set nchunks 16


poetry run -C ../../ceqtl -- \
    ceqtl med --name MediationAnalysis \
        --expr $exprfile \
        --geno $genofile \
        --trios $ceqtl_trios \
        --forks $nchunks \
        --nchunks $nchunks \
        --ComposeInputFiles.envs.transpose_expr \
        --ComposeInputFiles.envs.transpose_geno \
        --ComposeInputFiles.envs.transpose_cov \
        --cache force


poetry run -C ../../ceqtl -- \
    ceqtl med --name MediationAnalysis-cov \
        --expr $exprfile \
        --geno $genofile \
        --cov $covfile \
        --trios $ceqtl_trios_cov \
        --forks $nchunks \
        --nchunks $nchunks \
        --ComposeInputFiles.envs.transpose_expr \
        --ComposeInputFiles.envs.transpose_geno \
        --ComposeInputFiles.envs.transpose_cov \
        --profile slurm-med-lg \
        --cache force
