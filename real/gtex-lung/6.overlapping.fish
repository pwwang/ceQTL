#!/usr/bin/env fish

set ceqtl_genevar_pairs ./ceQTL-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set ceqtl_genevar_pairs_cov ./ceQTL-cov-output/CombinePValuesForVariantGenePairs/ceqtl-vargenes.txt
set goldfile ./eQTL-output/MatrixEQTL/cisqtls.txt
set goldfile_cov ./eQTL-cov-output/MatrixEQTL/cisqtls.txt

poetry run -C ../../ceqtl -- \
    python ./6.overlapping.py \
        --name eQTLceQTLVariantOverlapping \
        --eqtls $goldfile \
        --ceqtls $ceqtl_genevar_pairs \
        --eqtl-cols SNP \
        --ceqtl-cols SNP

poetry run -C ../../ceqtl -- \
    python ./6.overlapping.py \
        --name eQTLceQTLVariantOverlapping-cov \
        --eqtls $goldfile_cov \
        --ceqtls $ceqtl_genevar_pairs_cov \
        --eqtl-cols SNP \
        --ceqtl-cols SNP

poetry run -C ../../ceqtl -- \
    python ./6.overlapping.py \
        --name eQTLceQTLVgpairOverlapping \
        --eqtls $goldfile \
        --ceqtls $ceqtl_genevar_pairs

poetry run -C ../../ceqtl -- \
    python ./6.overlapping.py \
        --name eQTLceQTLVgpairOverlapping-cov \
        --eqtls $goldfile_cov \
        --ceqtls $ceqtl_genevar_pairs_cov
