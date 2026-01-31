#!/usr/bin/env fish

set ceqtl_genevar_pairs ./ceQTL-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set ceqtl_genevar_pairs_cov ./ceQTL-cov-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set goldfile ./eQTL-output/alleqtls.filtered.txt
set goldfile_cov ./eQTL-cov-output/alleqtls.filtered.txt

# set cutoff 1e-5
set cutoff 0.05

poetry run -C ../../ceqtl -- \
    python ../gtex-lung/6.overlapping.py \
        --name eQTLceQTLVariantOverlapping \
        --eqtls $goldfile \
        --ceqtls $ceqtl_genevar_pairs \
        --eqtl-cols SNP \
        --ceqtl-cols SNP \
        --eqtl-pos "FDR < $cutoff" \
        --ceqtl-pos "MetaPadj < $cutoff"

poetry run -C ../../ceqtl -- \
    python ../gtex-lung/6.overlapping.py \
        --name eQTLceQTLVariantOverlapping-cov \
        --eqtls $goldfile_cov \
        --ceqtls $ceqtl_genevar_pairs_cov \
        --eqtl-cols SNP \
        --ceqtl-cols SNP \
        --eqtl-pos "FDR < $cutoff" \
        --ceqtl-pos "MetaPadj < $cutoff"

poetry run -C ../../ceqtl -- \
    python ../gtex-lung/6.overlapping.py \
        --name eQTLceQTLVgpairOverlapping \
        --eqtls $goldfile \
        --ceqtls $ceqtl_genevar_pairs \
        --eqtl-pos "FDR < $cutoff" \
        --ceqtl-pos "MetaPadj < $cutoff"

poetry run -C ../../ceqtl -- \
    python ../gtex-lung/6.overlapping.py \
        --name eQTLceQTLVgpairOverlapping-cov \
        --eqtls $goldfile_cov \
        --ceqtls $ceqtl_genevar_pairs_cov \
        --eqtl-pos "FDR < $cutoff" \
        --ceqtl-pos "MetaPadj < $cutoff"
