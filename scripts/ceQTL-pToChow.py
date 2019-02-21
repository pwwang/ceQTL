"""
Prepare 3 files for chow-test.
1. Expression file
```
		X1  X2  X3  X4 ... Y
	S1  1   2   1   4  ... 9
	S2  2   3   1   1  ... 3
	... ...
	Sm  3   9   1   7  ... 8
```
2. Group file
```
        Case1	Case2
  S1	Group1	Group1
  S2	Group1	NA          # Exclude S2 in Case2
  S3	Group2	Group1
  ... ...
  Sm	Group2	Group2
```
3. Case file
```
  Case1	Y ~ X1
  Case2	Y ~ X2
```
"""
from collections import defaultdict
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

expfile   = {{i.expfile | quote}}
mutfile   = {{i.mutfile | quote}}
interfile = {{i.interfile | quote}}
covfile   = {{i.covfile | quote}}
genes     = {{i.genes.split(',') | repr}}
tfs       = {{i.tfs.split(';') | repr}}
outdata   = {{o.outdata | quote}}
outgroup  = {{o.outgroup | quote}}
outcase   = {{o.outcase | quote}}
genetfs   = { g:tfs[i].split(',') for i, g in enumerate(genes) }

# get gene, snp pairs
"""
chr1	12463073	12463074	AADACL4	0	+	chr1	12463073	12463074	chr1_12463073_rs7547740_A_G	0	+
chr1	12480504	12480505	AADACL4	0	+	chr1	12480504	12480505	chr1_12480504_rs6660365_T_C	0	+
chr1	12496021	12496022	AADACL4	0	+	chr1	12496021	12496022	chr1_12496021_rs6541023_T_C	0	+
"""
mutgenes   = defaultdict(lambda: [])
intereader = TsvReader(interfile)
genes = set()
for r in intereader:
	if not r[3] in genetfs:
		continue
	mutgenes[r[9]].append(r[3])
	genes.add(r[3])
intereader.close()

# shrink the sets
genetfs = {g: genetfs[g] for g in genes}
tfs     = list({tf for gtfs in genetfs.values() for tf in gtfs})

# nothing, write empty files
if not mutgenes or not genes:
	open(outdata, 'w').close()
	open(outgroup, 'w').close()
	open(outcase, 'w').close()
	exit(0)

# save the data file
# expfile
"""
	S1	S2	..	Sn
G1	...
G2	...
"""
expreader  = TsvReader(expfile)
expdata    = [r for r in expreader if r[0] in genes or r[0] in tfs]
expreader.close()
datawriter = TsvWriter(outdata)
for i, cname in enumerate(expreader.cnames):
	if i == 0:
		# genes + tfs
		datawriter.cnames = [r[0] for r in expdata]
		datawriter.writeHead()
	else:
		datawriter.write([cname] + [r[i] for r in expdata])
datawriter.close()
del expdata
genes = [g for g in genes if g in datawriter.cnames]
tfs   = [g for g in tfs if g in datawriter.cnames]

genetfs = {g: [tf for tf in gtfs if tf in tfs] for g, gtfs in genetfs.items() if g in genes}

# save the group file
# mutfile
"""
	S1	S2	..	Sn
M1	... (0/1/2/NA)
M2	...
"""
mutreader = TsvReader(mutfile)
mutdata   = [r for r in mutreader if r[0] in mutgenes]
mutreader.close()
groupwriter = TsvWriter(outgroup)
mutkids = lambda mut: {tf + '.' + g for g in mutgenes[mut] for tf in genetfs[g] if g in genes and tf != g}
mutrs   = lambda mut: mut.split('_')[2]
def mutgt(mut, gt):
	if gt == 'NA': return gt
	ref, alt = mut.split('_')[-2:]
	if gt == '2': return alt + alt + '(hom)'
	if gt == '1': return ref + alt + '(het)'
	return ref + ref + '(ref)'
for i, cname in enumerate(mutreader.cnames):
	if i == 0:
		groupwriter.cnames = ['.'.join([mutrs(r[0]), tfg]) for r in mutdata for tfg in mutkids(r[0])]
		groupwriter.writeHead()
	else:
		groupwriter.write([cname] + sum(([mutgt(r[0], r[i])] * len(mutkids(r[0])) for r in mutdata), []))
groupwriter.close()
del mutdata
del genetfs

# save the case file
# get the covariants
covariants = []
if (covfile):
	covreader  = TsvReader(covfile)
	covariants = covreader.cnames[1:]
	covreader.close()

bQuote = lambda s: '`{}`'.format(s)
casewriter = TsvWriter(outcase)
for case in groupwriter.cnames:
	(_, tf, g) = case.split('.')
	casewriter.write([case, "{} ~ {}".format(
		bQuote(g),
		' + '.join(bQuote(s) for s in ([tf] + covariants))
	)])
casewriter.close()

