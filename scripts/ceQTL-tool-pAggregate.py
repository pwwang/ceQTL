from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
method   = {{args.method | quote}}

def aggregate(pvals, method):
	if method == 'min':
		return min(pvals)
	elif method == 'fisher':
		from operator import mul
		try:
			reduce
		except NameError:
			from functools import reduce
		return reduce(mul, pvals, 1) ** (1.0/float(len(pvals)))
	else:
		raise ValueError('Method %s not supported yet.' % method)

def numpval(pval):
	try:
		return float(pval)
	except TypeError:
		return 1.0

reader    = TsvReader(infile)
writer    = TsvWriter(outfile)
prevsnp   = None
prevpvals = []
for r in reader:
	snp = r.Case.split('.')[0]
	if snp != prevsnp:
		if prevsnp:
			writer.write([
				prevsnp,
				aggregate(prevpvals, method)
			])
		prevsnp   = snp
		prevpvals = [numpval(r.Pval)]
	else:
		prevpvals.append(numpval(r.Pval))
writer.write([
	prevsnp,
	aggregate(prevpvals, method)
])

writer.close()
