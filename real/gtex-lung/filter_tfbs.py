#!/usr/bin/env python
import sys
from pathlib import Path

var_coord_file = ".pipen/ProcessGT/VariantCoords/0/output/variant_coords.bed"
tfbs_dir = "TFLINK_TFBS/chroms"

if len(sys.argv) < 2:
    print("Usage: filter_tfbs.py <input_file>")
    sys.exit(1)

print("Reading variant coordinates...")
var_coords = {}
with open(var_coord_file, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        snp = fields[3]
        if "rs" not in snp:
            continue

        var_coords[snp] = (fields[0], int(fields[2]))

infile = Path(sys.argv[1])
outfile = infile.with_suffix(".tfbs.txt")

with infile.open("r") as f:
    num_lines = sum(1 for line in f)

with infile.open("r") as fin, outfile.open("w") as fout:
    header = next(fin)
    fout.write(header)
    header_items = header.strip().split("\t")
    snp_index = header_items.index("SNP")
    tf_index = header_items.index("TF")
    i = 0
    for line in fin:
        if (i + 1) % 1000 == 0:
            print(f"Processed {i+1}/{num_lines} lines...")

        fields = line.strip().split("\t")
        snp = fields[snp_index]
        tf = fields[tf_index]
        if "rs" in snp:
            chrom, pos = var_coords[snp]
        else:
            chrom, pos = snp.split("_")[:2]
            pos = int(pos)

        tfbs_file = Path(tfbs_dir) / f"{tf}_{chrom}.bed"
        if not tfbs_file.exists():
            # print(f"Warning: {tfbs_file} not found.")
            i += 1
            continue

        with tfbs_file.open("r") as f:
            for tline in f:
                start, end = tline.strip().split("\t")[1:3]
                # print(start, end)
                start, end = int(float(start)), int(float(end))

                if start < pos <= end:
                    fout.write(line)
                    break

        i += 1
