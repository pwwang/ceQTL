#!/usr/bin/env fish

set exprfile (realpath .pipen/ProcessExpr/ConvertToGeneSymbol/0/output/Lung.v8.normalized_expression.matrix.txt)
set genofile (realpath ProcessGT-output/Plink2GTMat/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze-gtmat.txt)
set covfile (realpath ../gtex-lung/raw/Lung.v8.covariates.txt)

set program_name (basename (status --current-filename))
# if less than 4 arguments, print help
if test (count $argv) -lt 3
    echo "Usage: $program_name <SNP> <TF> <Target> [--cov]"
    exit 1
end

if contains -- --cov -- $argv
    set has_cov 1
else
    set has_cov 0
end

# remove --cov if any from argv, then get snp, tf, and target
set argv (string replace -- --cov '' $argv)
# remove empty strings
set argv (string split ' ' -- (string join ' ' -- $argv))
set snp $argv[1]
set tf $argv[2]
set target $argv[3]
set outdir (realpath "DensityPlot-output/Density-$snp-$tf-$target")

set cmd poetry run -C ../../ceqtl -- \
    ceqtl density \
        --snp $snp --tf $tf --target $target \
        --expr $exprfile --geno $genofile \
        --transpose-expr \
        --transpose-geno \
        --transpose-cov \
        --outdir $outdir

if test $has_cov -eq 1
    set cmd $cmd --cov $covfile
end

echo $cmd
command $cmd
