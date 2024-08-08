#!/usr/bin/env fish

set exprfile ProcessExpr-output/FilterGenes/Lung.v8.normalized_expression.matrix.txt
set genofile ProcessGT-output/Plink2GTMat/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.txt
set snppos ProcessGT-output/VariantCoords/variant_coords.bed
set genepos ~/reference/hg38/hg38.refGene.transcript.gtf
set covfile raw/Lung.v8.covariates.txt
set dist 2500


poetry run -C ../../ceqtl -- \
    ceqtl eqtl --name eQTL \
    --MatrixEQTL.in.geno $genofile --MatrixEQTL.in.expr $exprfile \
    --MatrixEQTL.envs.snppos $snppos --MatrixEQTL.envs.genepos $genepos \
    --MatrixEQTL.envs.dist $dist

poetry run -C ../../ceqtl -- \
    ceqtl eqtl --name eQTL-cov \
    --MatrixEQTL.in.geno $genofile --MatrixEQTL.in.expr $exprfile --MatrixEQTL.in.cov $covfile \
    --MatrixEQTL.envs.snppos $snppos --MatrixEQTL.envs.genepos $genepos \
    --MatrixEQTL.envs.dist $dist
