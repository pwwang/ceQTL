#!/usr/bin/env fish

set exprfile ProcessExpr-output/FilterGenes/Lung.v8.normalized_expression.matrix.txt
set genofile ProcessGT-output/Plink2GTMat/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.txt
set tftarget raw/TF-Gene.gmt
set genesnp raw/gene_var.gmt
set triofile raw/TF-snp-TG.txt
set covfile ../gtex-lung/raw/Lung.v8.covariates.txt
set forks 16
set nchunks 1
set ncores 16
set config ./3.call.toml

poetry run -C ../../ceqtl -- \
	ceqtl call @$config --name ceQTL \
		--forks $forks --nchunks $nchunks --ncores $ncores \
		--expr $exprfile --geno $genofile \
		--tftarget $tftarget --genesnp $genesnp \
		--triofile $triofile \
		--DataPreparation.envs.transpose_geno \
		--DataPreparation.envs.transpose_expr

poetry run -C ../../ceqtl -- \
	ceqtl call @$config --name ceQTL-cov \
		--forks $forks --nchunks $nchunks --ncores $ncores \
		--expr $exprfile --geno $genofile \
		--tftarget $tftarget --genesnp $genesnp --cov $covfile \
		--triofile $triofile \
		--DataPreparation.envs.transpose_cov \
		--DataPreparation.envs.transpose_geno \
		--DataPreparation.envs.transpose_expr
