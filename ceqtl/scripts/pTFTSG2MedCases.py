
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

tftfile = {{i.tftfile | quote}}
sgfile = {{i.sgfile | quote}}
outfile = {{o.outfile | quote}}

# tftfile:
#  TF  GENE

# sgfile:
#  SNP GENE

# outfile:
#  CaseX    GENE ~ TF | SNP

tft_reader = TsvReader(tftfile, cnames=False)
tft = {}
for row in tft_reader:
    tft.setdefault(row[1], []).append(row[0])
tft_reader.close()

sg_reader = TsvReader(sgfile, cnames=False)
sg = {}
for row in sg_reader:
    sg.setdefault(row[1], []).append(row[0])
sg_reader.close()

writer = TsvWriter(outfile)
index = 1
for gene in set(tft) & set(sg):
    for tf in tft[gene]:
        for snp in sg[gene]:
            writer.write([
                'Case%s' % index,
                '`%s` ~ `%s` | `%s`' % (gene, tf, snp)
            ])
            index += 1
writer.close()
