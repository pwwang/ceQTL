#!/usr/bin/env python
from __future__ import print_function
from os import path
from ceQTL import VERSION
from pyppl import Proc, PyPPL
from pyppl.utils import split
from pyppl.parameters import commands
from bioprocs import params
if path.isfile('./config.yaml'):
	params.loadFile('./config.yaml')

from bioprocs.common import pSort, pFile2Proc, pStr2File
from bioprocs.snp import pRs2Bed
from bioprocs.tsv import pTsvJoin, pTsv, pTsvMerge
from bioprocs.bed import pBedSort, pBedIntersect
from bioprocs.plot import pManhattan
from bioprocs.tfbs import pAtSnp
from bioprocs.vcfnext import pGTMat2Bed
from bioprocs.utils import shell

commands._desc = [
	'ceQTL-tool.py v' + VERSION,
	'A set of tools for ceQTL.py'
]
commands.version        = 'Show current version of ceQTL.py'
commands.version._hbald = False

commands.manhattan           = 'Do manhattan plot for ceQTLs'
commands.man                 = commands.manhattan
commands.man.cefile.desc     = 'The output from ceQTL pipeline.'
commands.man.cefile.required = True
commands.man.aggr            = 'min'
commands.man.aggr.desc       = 'Method to aggregate pvalues (min or fisher).'
commands.man.outfile         = './ceqtl.manhattan.png'
commands.man.outfile.desc    = 'The output file for the manhattan plot.'
commands.man.hifile          = ''
commands.man.hifile.desc     = [
	'A file with list of snps you want to highlight in manhattan plot.', 'One per line.'
]

commands.atsnp              = 'Calculate the differential binding affinitiy w/o the mutation.'
commands.atsnp.cefile       = commands.man.cefile
commands.atsnp.nthread      = 1
commands.atsnp.nthread.desc = 'Number threads to use.'
commands.atsnp.tflist       = params.tflist
commands.atsnp.tflist.show  = True
commands.atsnp.motifdb      = params.tfmotifs
commands.atsnp.motifdb.show = True
commands.atsnp.outfile      = './ceqtl.atsnp.txt'
commands.atsnp.outfile.desc = 'The output file of atSNP.'
commands.atsnp.man.desc     = 'File of the manhattan plot for the SNPs, `False/None` to disable.'
commands.atsnp.hifile       = commands.man.hifile

commands.filter            = 'Filter the ceqtl results.'
commands.filter.cefile     = commands.man.cefile
commands.filter.snps.desc  = 'Filter by a list of SNPs, or a file with SNPs, one per line.'
commands.filter.snps.type  = list
commands.filter.genes.desc = 'Filter by a list of genes, or a file with genes, one per line.'
commands.filter.genes.type = list
commands.filter.tfs.desc   = 'Filter by a list of tfs, or a file with tfs, one per line.'
commands.filter.tfs.type   = list
commands.filter.regs.desc  = 'Filter by a list of regions (in format of chrA:POS1-POS2) or a BED file.'
commands.filter.regs.type  = list
commands.filter.row.desc   = [
	'A python lambda function with row to apply filters on.',
	'For example, "lambda row: float(row.Pval) < 0.05" will keep records with Pval < 0.05',
	'Available columns: [Case, Pooled, Groups, Fstat, Pval] and maybe "AdjPval",',
	'If the ceQTLs are generated with "-padj" is on.'
]
commands.filter.helper       = ''
commands.filter.helper.desc  = 'Helper code for row filter.'
commands.filter.outfile      = './ceqtl.filtered.txt'
commands.filter.desc         = 'The output filtered ceqtl file.'
commands.filter.connect      = 'and'
commands.filter.connect.desc = 'How are those filtered connected (and/or)'
commands.filter.row.callback = lambda opt, args: None if args.snps.value or args.genes.value or args.tfs.value or args.regs.value or opt.value else 'One of filters is required: snps, genes, tfs, regs or row.'

### common processes	
# Case	Pooled	Groups	Fstat	Pval	AdjPval
# chr12_125754313_rs10744198_C_G.ZNF219.SPPL3	ZNF219=-0.257,_=8.809,N=703	AA(hom):ZNF219=-0.163,_=8.422,N=115; GA(het):ZNF219=-0.404,_=9.406,N=332; GG(ref):ZNF219=-0.078,_=8.089,N=256	4.754	0.000865	9.977E-01
pSortByRs                     = pSort.copy()
pSortByRs.args.inopts.skip    = 1
pSortByRs.args.inopts.delimit = '_'
pSortByRs.args.params.k       = 3

# including coordinates
pSortBySnp                     = pSort.copy()
pSortBySnp.args.inopts.skip    = 1

pSortByTF                     = pSort.copy()
pSortByTF.args.inopts.skip    = 1
pSortByTF.args.inopts.delimit = '.'
pSortByTF.args.params.k       = 2

pSortByGene                     = pSort.copy()
pSortByGene.args.inopts.skip    = 1
pSortByGene.args.inopts.delimit = '.'
pSortByGene.args.params.k       = 3

def setOutfile(p, outfile, key = 'outfile'):
	p.exdir  = path.dirname(outfile)
	p.expart = key
	outputs  = p.config['output']
	if isinstance(outputs, list):
		p.config['output'][0] = '{}:file:{}'.format(key, path.basename(outfile))
	else:
		outputs = split(outputs, ',')
		outputs[0] = '{}:file:{}'.format(key, path.basename(outfile))
		p.config['output'] = outputs

def ceQTL_manhattan(args):

	pSortByRs.input               = [args.cefile]

	pAggregate               = Proc(desc = 'Aggregate pvalues each snp.')
	pAggregate.depends       = pSortByRs
	pAggregate.input         = 'infile:file'
	pAggregate.output        = 'outfile:file:{{i.infile | bn}}'
	pAggregate.args.method   = args.aggr
	pAggregate.lang          = params.python.value
	pAggregate.script        = 'file:scripts/ceQTL-tool-pAggregate.py'

	pGTMat2Bed.depends = pAggregate
	pGTMat2Bed.args.inopts.cnames = False
	pGTMat2Bed.args.ncol = 66

	pBedSort.depends       = pGTMat2Bed
	pBedSort.args.chrorder = params.chrorder.value

	pManhattan.depends = pBedSort
	if args.hifile:
		pManhattan.input = lambda ch: ch.cbind(args.hifile)
	pManhattan.args.gsize = params.gsize.value
	setOutfile(pManhattan, args.outfile)

	PyPPL().start(pSortByRs).run()

def ceQTL_atsnp(args):
	
	pGTMat2Bed.input               = [args.cefile]
	pGTMat2Bed.args.inopts.cnames  = False
	pGTMat2Bed.args.inopts.skip    = 1
	pGTMat2Bed.args.inopts.delimit = '.'
	pGTMat2Bed.args.name           = 'full'
	pGTMat2Bed.args.ncol           = 8

	pSortSnp             = pSort.copy()
	pSortSnp.depends     = pGTMat2Bed
	pSortSnp.args.unique = True

	pSortByTF.input = [args.cefile]

	pFilterTFs                     = pTsvJoin.copy()
	pFilterTFs.depends             = pSortByTF
	pFilterTFs.input               = lambda ch: [ch.insert(0, args.tflist).flatten()]
	pFilterTFs.args.inopts.cnames  = False
	pFilterTFs.args.inopts.skip    = [0, 1]
	pFilterTFs.args.inopts.delimit = ['\t', '.']
	pFilterTFs.args.outopts.cnames = False
	pFilterTFs.args.match          = 'lambda r1, r2: TsvJoin.compare(r1[1], r2[1])'
	pFilterTFs.args.do             = 'lambda out, r1, r2: out.write(r1)'

	pAtSnp.depends       = pFilterTFs, pSortSnp
	pAtSnp.args.tfmotifs = args.motifdb
	pAtSnp.args.fdr      = False
	pAtSnp.args.plot     = False
	pAtSnp.args.nthread  = args.nthread
	setOutfile(pAtSnp, args.outfile)

	if args.man:
		pToMan                     = pTsv.copy()
		pToMan.depends             = pAtSnp
		pToMan.args.outopts.cnames = False
		# [chr1, 12496021, rs6541023, 0.04]
		pToMan.args.helper         = 'snprec = lambda x: [x[0], int(x[1]) - 1, x[1], x[2], 0, "+"]'
		pToMan.args.row            = 'lambda r: snprec(r.Snp.split("_")[:3]) + [r.Pval_Diff]'
		
		pBedSort.depends       = pToMan
		pBedSort.args.chrorder = params.chrorder.value

		pManhattan.depends = pBedSort
		if args.hifile:
			pManhattan.input = lambda ch: ch.cbind(args.hifile)
		pManhattan.args.gsize = params.gsize.value
		setOutfile(pManhattan, args.man)

	PyPPL().start(pGTMat2Bed, pSortByTF).run()

def ceQTL_filter_snps(args, pCefile, hasAdjPval, starts):
	if not args.snps:
		return pCefile
	header = ['Case', 'Pooled', 'Groups', 'Fstat', 'Pval']
	if hasAdjPval:
		header.append('AdjPval')

	if not path.isfile(args.snps[0]):
		pSnpfile       = pStr2File.copy()
		pSnpfile.input = [','.join(sorted(args.snps))]
	else:
		pSnpfile = pFile2Proc.copy()
		pSnpfile.input = [args.snps[0]]
	starts.append(pSnpfile)
	pSortByRs.depends = pCefile

	pFilterSnps                           = pTsvJoin.copy()
	pFilterSnps.depends                   = pSortByRs,  pSnpfile
	pFilterSnps.input                     = lambda ch1, ch2: [[ch1.get(), ch2.get()]]
	pFilterSnps.args.inopts.cnames        = False
	pFilterSnps.args.inopts.skip          = [1, 0]
	pFilterSnps.args.inopts.delimit       = ['_', '\t']
	pFilterSnps.args.outopts.cnames       = header
	pFilterSnps.args.outopts.delimit      = '_'
	pFilterSnps.args.outopts.headCallback = 'lambda cnames: "\t".join(cnames)'
	pFilterSnps.args.match                = 'lambda r1, r2: TsvJoin.compare(r1[2], r2[0])'
	pFilterSnps.args.do                   = 'lambda out, r1, r2: out.write(r1.values())'

	return pFilterSnps

def ceQTL_filter_genes(args, pCefile, hasAdjPval, starts):
	if not args.genes:
		return pCefile
	header = ['Case', 'Pooled', 'Groups', 'Fstat', 'Pval']
	if hasAdjPval:
		header.append('AdjPval')

	if not path.isfile(args.genes[0]):
		pGenefile       = pStr2File.copy()
		pGenefile.input = [','.join(sorted(args.genes))]
	else:
		pGenefile = pFile2Proc.copy()
		pGenefile.input = [args.genes[0]]
	starts.append(pGenefile)
	pSortByGene.depends = pCefile

	pFilterGenes                           = pTsvJoin.copy()
	pFilterGenes.depends                   = pSortByGene, pGenefile
	pFilterGenes.input                     = lambda ch1, ch2: [[ch1.get(), ch2.get()]]
	pFilterGenes.args.inopts.cnames        = False
	pFilterGenes.args.inopts.skip          = [1, 0]
	pFilterGenes.args.inopts.delimit       = ['.', '\t']
	pFilterGenes.args.outopts.cnames       = header
	pFilterGenes.args.outopts.delimit      = '.'
	pFilterGenes.args.outopts.headCallback = 'lambda cnames: "\\t".join(cnames)'
	pFilterGenes.args.match                = 'lambda r1, r2: TsvJoin.compare(r1[2].split("\\t")[0], r2[0])'
	pFilterGenes.args.do                   = 'lambda out, r1, r2: out.write(r1.values())'

	return pFilterGenes

def ceQTL_filter_tfs(args, pCefile, hasAdjPval, starts):
	if not args.tfs:
		return pCefile
	header = ['Case', 'Pooled', 'Groups', 'Fstat', 'Pval']
	if hasAdjPval:
		header.append('AdjPval')

	if not path.isfile(args.tfs[0]):
		pTffile       = pStr2File.copy()
		pTffile.input = [','.join(sorted(args.tfs))]
	else:
		pTffile = pFile2Proc.copy()
		pTffile.input = [args.tfs[0]]
	starts.append(pTffile)
	pSortByTF.depends = pCefile

	pFilterTfs                           = pTsvJoin.copy()
	pFilterTfs.depends                   = pSortByTF,  pTffile
	pFilterTfs.input                     = lambda ch1, ch2: [[ch1.get(), ch2.get()]]
	pFilterTfs.args.inopts.cnames        = False
	pFilterTfs.args.inopts.skip          = [1, 0]
	pFilterTfs.args.inopts.delimit       = ['.', '\t']
	pFilterTfs.args.outopts.cnames       = header
	pFilterTfs.args.outopts.delimit      = '.'
	pFilterTfs.args.outopts.headCallback = 'lambda cnames: "\t".join(cnames)'
	pFilterTfs.args.match                = 'lambda r1, r2: TsvJoin.compare(r1[1], r2[0])'
	pFilterTfs.args.do                   = 'lambda out, r1, r2: out.write(r1.values())'

	return pFilterTfs

def ceQTL_filter_regs(args, pCefile, hasAdjPval, starts):
	if not args.regs:
		return pCefile
	header = ['Case', 'Pooled', 'Groups', 'Fstat', 'Pval']
	if hasAdjPval:
		header.append('AdjPval')

	if not path.isfile(args.regs[0]):
		pRegfile       = pStr2File.copy()
		pRegfile.input = [','.join(reg.replace(':', '\t').replace('-', '\t') for reg in args.regs)]
	else:
		pRegfile = pFile2Proc.copy()
		pRegfile.input = [args.regs[0]]
	starts.append(pRegfile)

	pGTMat2Bed.depends             = pCefile
	pGTMat2Bed.args.inopts.cnames  = False
	pGTMat2Bed.args.inopts.skip    = 1
	pGTMat2Bed.args.inopts.delimit = '.'
	pGTMat2Bed.args.name           = 'full'
	pGTMat2Bed.args.ncol           = 6

	pUniqueBed             = pSort.copy()
	pUniqueBed.depends     = pGTMat2Bed
	pUniqueBed.args.unique = True

	pBedIntersect.depends = pUniqueBed, pRegfile
	pSortBySnp.depends = pCefile

	pFilterRegs                           = pTsvJoin.copy()
	pFilterRegs.depends                   = pSortBySnp,  pBedIntersect
	pFilterRegs.input                     = lambda ch1, ch2: [[ch1.get(), ch2.get()]]
	pFilterRegs.args.inopts.cnames        = False
	pFilterRegs.args.inopts.skip          = [1, 0]
	pFilterRegs.args.inopts.delimit       = ['.', '\t']
	pFilterRegs.args.outopts.cnames       = header
	pFilterRegs.args.outopts.delimit      = '.'
	pFilterRegs.args.outopts.headCallback = 'lambda cnames: "\t".join(cnames)'
	pFilterRegs.args.match                = 'lambda r1, r2: TsvJoin.compare(r1[0], r2[3])'
	pFilterRegs.args.do                   = 'lambda out, r1, r2: out.write(r1.values())'
	
	return pFilterRegs
	

def ceQTL_filter_row(args, pCefile, hasAdjPval, starts):
	if not args.row:
		return pCefile
	
	pFilterRow = pTsv.copy()
	pFilterRow.depends = pCefile
	pFilterRow.args.helper = args.helper
	pFilterRow.args.row = args.row
	if hasAdjPval:
		pFilterRow.args.inopts.row = 'lambda r: setattr(r, "Pval", float(r.Pval)) or setattr(r, "AdjPval", float(r.AdjPval)) or r'
	else:
		pFilterRow.args.inopts.row = 'lambda r: setattr(r, "Pval", float(r.Pval)) or r'
	return pFilterRow

def ceQTL_filter(args):
	# check if AdjPval is in the header
	p = shell.Shell().head(n = 1, _ = args.cefile).pipe().grep('AdjPval').run(raiseExc = False, logger = False)

	hasAdjPval    = p.rc == 0
	pCefile       = pFile2Proc.copy()
	pCefile.input = [args.cefile]
	starts        = [pCefile]

	pCefile_snps  = ceQTL_filter_snps (args, pCefile, hasAdjPval, starts)
	pCefile_genes = ceQTL_filter_genes(args, (pCefile, pCefile_snps)[int(args.connect == 'and')], hasAdjPval, starts)
	pCefile_tfs   = ceQTL_filter_tfs  (args, (pCefile, pCefile_genes)[int(args.connect == 'and')], hasAdjPval, starts)
	pCefile_regs  = ceQTL_filter_regs (args, (pCefile, pCefile_tfs)[int(args.connect == 'and')], hasAdjPval, starts)
	pCefile_row   = ceQTL_filter_row  (args, (pCefile, pCefile_regs)[int(args.connect == 'and')], hasAdjPval, starts)

	if args.connect == 'and':
		setOutfile(pCefile_row, args.outfile)
	else:
		procs = list({p for p in {pCefile_snps, pCefile_genes, pCefile_tfs, pCefile_regs, pCefile_row} if not p is pCefile})
		pTsvMerge.depends = procs
		pTsvMerge.input   = lambda *chs: [sum((ch.flatten() for ch in chs), [])]
		pTsvMerge.args.inopts.cnames = True
		
		pSortMerged             = pSort.copy()
		pSortMerged.depends     = pTsvMerge
		pSortMerged.args.unique = True

		setOutfile(pSortMerged, args.outfile)

	PyPPL().start(starts).run()

def ceQTL_version():
	print('ceQTL-tool.py v' + VERSION)

if __name__ == "__main__":
	command, args = commands.parse()
	if command == 'version':
		ceQTL_version()
	elif command in ['man', 'manhattan']:
		ceQTL_manhattan(args)
	elif command == 'atsnp':
		ceQTL_atsnp(args)
	elif command == 'filter':
		ceQTL_filter(args)
