"""
Convert Snp-gene pairs to bed files to decide cis-eQTL calculation
"""
import math
import random
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

# snp gene
# SNP1 Gene10 # sorted by gene
infile = {{i.infile | quote}}
snpfile = {{o.snpfile | quote}}
genefile = {{o.genefile | quote}}
snppergene = {{args.snppergene | repr}}
nchr = {{args.nchr | repr}}
seed = {{args.seed | repr}}
# distances between genes
dist = {{args.dist | repr}}

random.seed(seed)

reader = TsvReader(infile, cnames=False)
allsnps = set(reader.dump(0))
reader.rewind()
allgenes = set(reader.dump(1))
reader.close()

# assign a probability to each snp
nsnps = len(allsnps)
ngenes = len(allgenes)
snp_probs = dict(zip(allsnps, random.choices(range(ngenes * snppergene),
                                             k=nsnps)))

genebed = TsvWriter(genefile)
snpbed = TsvWriter(snpfile)

geneperchr = math.ceil(float(ngenes) / float(nchr))
for i, gene in enumerate(allgenes):
    chrname = 'chr' + str(int(i % nchr) + 1)
    start = (int(i / nchr) + 1) * dist
    end = start + 1
    first_snp_pos = int(start - dist/2.0 - snppergene)
    snps = (snp for snp in snp_probs
            if i * snppergene <= snp_probs[snp] < (i+1)*snppergene)
    genebed.write([chrname, start, end, gene, 0, '+'])
    for j, snp in enumerate(snps):
        snppos = first_snp_pos + j
        snpbed.write([chrname, snppos, snppos, snp, 0, '+'])

genebed.close()
snpbed.close()
