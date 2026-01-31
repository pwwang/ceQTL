#!/usr/bin/env fish

set motiffile ./ceQTL-output/CombineAndAdjustPValues/ceqtl-trios.txt
set varfile ./ProcessGT-output/VariantCoords/variant_coords.bed
set motiffile_cov ./ceQTL-cov-output/CombineAndAdjustPValues/ceqtl-trios.sig.txt
set ncores 4

# poetry run -C ../../ceqtl -- \
#     pipen run regulatory MotifAffinityTest \
#     --name MotifAffinityTest \
#     --in.motiffile $motiffile \
#     --in.varfile $varfile \
#     --envs.ncores $ncores \
#     --envs.regulator_col TF \
#     --envs.tool atsnp \
#     --envs.notfound ignore \
#     --envs.genome hg38

poetry run -C ../../ceqtl -- \
    pipen run regulatory MotifAffinityTest \
    --name MotifAffinityTest-cov \
    --in.motiffile $motiffile_cov \
    --in.varfile (cat $varfile | grep -f (cut -f1 $motiffile_cov | psub) | psub) \
    --envs.ncores $ncores \
    --envs.regulator_col TF \
    --envs.tool atsnp \
    --envs.notfound ignore \
    --envs.genome hg38
