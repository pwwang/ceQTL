#!/usr/bin/env fish

set varfile ./ProcessGT-output/VariantCoords/variant_coords.bed
set motiffile ./ceQTL-output/CombineAndAdjustPValues/ceqtl-trios.txt
set motiffile_cov ./ceQTL-cov-output/CombineAndAdjustPValues/ceqtl-trios.txt
set ncores 4

if test (count $argv) -lt 2
    set program_name (basename (status --current-filename))
    echo "Usage: $program_name <SNP> <TF> [--cov]"
    exit 1
end

if contains -- --cov -- $argv
    set has_cov 1
else
    set has_cov 0
end

set snp $argv[1]
set tf $argv[2]


if test $has_cov -eq 1
    set outdir ./MotifAffinityTestIndividually/$snp-$tf-cov
    mkdir -p $outdir

    set tmp_motiffile $TMPDIR/7.affinity_test_individually/$snp-$tf-cov.txt
    mkdir -p (dirname $tmp_motiffile)
    head -1 $motiffile_cov > $tmp_motiffile
    grep "\s$snp\s$tf\s" $motiffile_cov | uniq >> $tmp_motiffile

    set tmp_varfile $TMPDIR/7.affinity_test_individually/$snp-$tf-cov.bed
    mkdir -p (dirname $tmp_varfile)
    grep "\s$snp\s" $varfile > $tmp_varfile

    set cmd poetry run -C ../../ceqtl -- \
        pipen run regulatory MotifAffinityTest \
        --name MotifAffinityTest_{$snp}_{$tf}-cov \
        --outdir $outdir \
        --in.motiffile $tmp_motiffile \
        --in.varfile $tmp_varfile \
        --envs.ncores $ncores \
        --envs.regulator_col TF \
        --envs.tool atsnp \
        --envs.notfound ignore \
        --envs.cutoff 1 \
        --envs.genome hg38
else
    set outdir ./MotifAffinityTestIndividually/$snp-$tf
    mkdir -p $outdir

    set tmp_motiffile $TMPDIR/7.affinity_test_individually/$snp-$tf.txt
    mkdir -p (dirname $tmp_motiffile)
    head -1 $motiffile > $tmp_motiffile
    grep "\s$snp\s$tf\s" $motiffile | uniq >> $tmp_motiffile

    set tmp_varfile $TMPDIR/7.affinity_test_individually/$snp-$tf.bed
    mkdir -p (dirname $tmp_varfile)
    grep "\s$snp\s" $varfile > $tmp_varfile

    set cmd poetry run -C ../../ceqtl -- \
        pipen run regulatory MotifAffinityTest \
        --name MotifAffinityTest_{$snp}_{$tf} \
        --outdir $outdir \
        --in.motiffile $tmp_motiffile \
        --in.varfile $tmp_varfile \
        --envs.ncores $ncores \
        --envs.regulator_col TF \
        --envs.tool atsnp \
        --envs.notfound ignore \
        --envs.cutoff 1 \
        --envs.genome hg38
end

echo $cmd
command $cmd
