#!/usr/bin/env fish

set exprfile .pipen/ceQTL/DataPreparation/0/output/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.prepared/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.chunk-1/expression-nocov.txt
set genofile .pipen/ceQTL/DataPreparation/0/output/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.prepared/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.chunk-1/genotype.txt
set snppos ProcessGT-output/VariantCoords/variant_coords.bed
set genepos ~/reference/hg38/hg38.refGene.transcript.gtf
set covfile ../gtex-lung/raw/Lung.v8.covariates.txt
# set dist 2500


poetry run -C ../../ceqtl -- \
    ceqtl eqtl --name eQTL \
    --MatrixEQTL.in.geno $genofile --MatrixEQTL.in.expr $exprfile \
    --MatrixEQTL.envs.transpose_geno \
    --MatrixEQTL.envs.transpose_expr \
    --MatrixEQTL.envs.transp 1.1

poetry run -C ../../ceqtl -- \
    ceqtl eqtl --name eQTL-cov \
    --MatrixEQTL.in.geno $genofile --MatrixEQTL.in.expr $exprfile --MatrixEQTL.in.cov $covfile \
    --MatrixEQTL.envs.transpose_geno \
    --MatrixEQTL.envs.transpose_expr \
    --MatrixEQTL.envs.transp 1.1

# filter the eQTLs with snp-gene pairs
set sgfile raw/snp-gene.txt
# ABCC1   chr10_100006504_T_C_b38 chr10_100267608_T_C_b38 chr10_100267650_C_T_b38
# ABCC3   chr10_100006504_T_C_b38 chr10_100267608_T_C_b38 chr10_100267650_C_T_b38
set eqtlfile eQTL-output/MatrixEQTL/alleqtls.txt
set eqtlfile_cov eQTL-cov-output/MatrixEQTL/alleqtls.txt
# SNP     gene  beta    t-stat  p-value FDR
# chr1_63523233_G_A_b38   ITGB3BP
# chr22_42183514_T_C_b38  CYP2D6
# chr22_41951056_T_C_b38  CYP2D6

set eqtlfile_filtered eQTL-output/alleqtls.filtered.txt
set eqtlfile_cov_filtered eQTL-cov-output/alleqtls.filtered.txt

head -1 $eqtlfile > $eqtlfile_filtered
grep -wf $sgfile $eqtlfile | sort -u >> $eqtlfile_filtered

head -1 $eqtlfile_cov > $eqtlfile_cov_filtered
grep -wf $sgfile $eqtlfile_cov | sort -u >> $eqtlfile_cov_filtered
