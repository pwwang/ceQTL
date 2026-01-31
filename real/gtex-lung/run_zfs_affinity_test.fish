#!/usr/bin/env fish

set motiffile (realpath ./zhifu.trio.091725.txt)
set varfile ./ProcessGT-output/VariantCoords/variant_coords.bed
set ncores 4

poetry run -C ../../ceqtl -- \
    pipen run regulatory MotifAffinityTest \
    --name MotifAffinityTest-zfs \
    --in.motiffile $motiffile \
    --in.varfile (cat $varfile | grep -f (cut -f3 $motiffile | psub) | psub) \
    --envs.ncores $ncores \
    --envs.regulator_col TF.ceQTL \
    --envs.var_col snp.ceQTL \
    --envs.tool atsnp \
    --envs.cutoff 1.1 \
    --envs.plot_nvars 100 \
    --envs.notfound ignore \
    --envs.atsnp_args '{"padj_cutoff": "FALSE"}' \
    --envs.genome hg38