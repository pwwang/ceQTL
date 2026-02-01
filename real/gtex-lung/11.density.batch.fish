#!/usr/bin/env fish

# Run SNP-TF-Target density plot from zfs_trio.txt
# SNP is at column: snp.eQTL
# TF is at column: TF.eQTLX
# Target is at column: TG.eQTLX

while read line
    set snp (echo $line | cut -d(echo -e "\t") -f4)
    set tf (echo $line | cut -d(echo -e "\t") -f9)
    set target (echo $line | cut -d(echo -e "\t") -f8)
    # echo "Processing SNP: $snp, TF: $tf, Target: $target"
    echo "fish 11.density.fish $snp $tf $target 2>&1 >/tmp/$snp-$tf-$target.density.batch.log"
end < zfs_trio.txt

# then run with GNU parallel or sequentially
# e.g.
# fish 11.density.batch.fish | parallel -j 20 --progress
