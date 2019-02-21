#!/usr/bin/env python
from os import path
from collections import defaultdict
from pyppl import PyPPL, Channel, logger, Proc
from bioprocs import params

if path.isfile('./config.yaml'):
	params.loadFile('./config.yaml')

params.runner          = 'sge1d'
params.runner.desc     = 'The runner for some heavy procs.'
params.forks           = 200
params.forks.desc      = 'How many batches to run simutanuously.'
params.tfhits.required = True
params.tfhits.desc     = [
	'The TF-gene pairs indicating the potential regulation relationship.', 
	'It could be a BED file with the name column "<TF>::<Gene>", or', 
	'a tab-delimited file with first column TFs and second column Genes.']
params.msmat.required  = True
params.msmat.desc      = 'The genotype matrix.'
params.expmat.required = True
params.expmat.desc     = 'The gene expression matrix.'
params.covfile         = ''
params.covfile.desc    = 'The covariate file.'
params.padj            = False
params.padj.desc       = [
	'Do adjustment on pvalues.',
	'  * `True`: Enable the action, and use method `fdr`',
	'  * `False`: Disable the action',
	'  * See `?p.adjust` in R for other available methods.']
params.dist            = 250000
params.dist.desc       = 'The up/down-stream distance to select variants'
params.pval            = 1e-3
params.pval.desc       = 'The p-value cutoff'
params.outfile         = './ceqtl.txt'
params.outfile.desc    = 'The output file'
params.ppldir          = './workdir'
params.ppldir.desc     = 'The workdir direcotry'
params.batch           = 1000
params.batch.desc      = 'Divide all genes into N batches for parallel computing'

from bioprocs.common import pSort
from bioprocs.gene import pPromoters
from bioprocs.bed import pBedIntersect
from bioprocs.tsv import pTsvMerge
from bioprocs.stats import pChow
from bioprocs.vcfnext import pGTMat2Bed
from bioprocs.utils.parallel import distributeList
from bioprocs.utils.tsvio2 import TsvReader

params = params.parse()

logger.logger.info('Reading tfhits ...')

# get tf-gene pairs
genes    = defaultdict(lambda: set())
reader   = TsvReader(params.tfhits)
allgenes = set()
npairs   = 0
if params.tfhits.endswith('.bed'):
	for r in reader:
		tf, gene = r[3].split('::')
		genes[gene].add(tf)
		allgenes.add(gene)
		npairs += 1
else:
	for r in reader:
		tf, gene = r[:2]
		genes[gene].add(tf)
		allgenes.add(gene)
		npairs += 1
reader.close()
logger.logger.info('Got %s genes and %s pairs.', len(allgenes), npairs)

splits = distributeList(list(allgenes), params.batch)
del allgenes

pExpmat                  = pSort.copy()
pExpmat.desc             = 'Sort expression matrix by gene names.'
pExpmat.input            = [params.expmat]
pExpmat.args.inopts.skip = 1

pMsmat                  = pSort.copy()
pMsmat.desc             = 'Sort mutation matrix by coordinates.'
pMsmat.input            = [params.msmat]
pMsmat.args.params.k    = ['1,1', '2,2n']
pMsmat.args.inopts.skip = 1

pGTMat2Bed.depends   = pMsmat
pGTMat2Bed.desc      = 'Convert mutation matrix to BED file'
pGTMat2Bed.args.name = 'full'
pGTMat2Bed.args.ncol = 6

pPromoters.depends              = pExpmat
pPromoters.args.inopts.cnames   = True
pPromoters.args.region.up       = params.dist
pPromoters.args.region.withbody = True

pBedIntersect.desc            = 'Overlap candicate regions with mutation.'
pBedIntersect.depends         = pPromoters, pGTMat2Bed
pBedIntersect.args.params.wb  = True
pBedIntersect.args.params.wao = False

pSortInter               = pSort.copy()
pSortInter.desc          = 'Sort mutation and gene intersect file.'
pSortInter.depends       = pBedIntersect
pSortInter.args.params.k = ['1,1', '2,2n']

pToChow         = Proc(desc = 'Prepare files for Chow test')
pToChow.input   = 'expfile:file, mutfile:file, interfile:file, covfile:file, genes, tfs'
pToChow.depends = pExpmat, pMsmat, pSortInter
pToChow.runner  = params.runner
pToChow.input   = lambda ch1, ch2, ch3: ch1.cbind(ch2, ch3).cbind(params.covfile).cbind(Channel([
	(
		','.join(sgenes), 
		';'.join(','.join(tf for tf in genes[sg]) for sg in sgenes)
	) for sgenes in splits
]))
pToChow.output = [
	'outdata:file:job{{job.index + 1}}.chowdata.txt',
	'outgroup:file:job{{job.index + 1}}.chowgroup.txt',
	'outcase:file:job{{job.index + 1}}.chowcase.txt',
]
pToChow.lang   = params.python
pToChow.script = 'file:scripts/ceQTL-pToChow.py'

pChow.depends   = pToChow
pChow.runner    = params.runner
pChow.args.plot = False
pChow.args.fdr  = False
pChow.args.cov  = params.covfile
pChow.args.pval = 1 if params.padj else params.pval 

pTsvMerge.depends            = pChow
pTsvMerge.input              = lambda ch: [ch.flatten(0)]
pTsvMerge.args.inopts.cnames = True
if not params.padj:
	pTsvMerge.exdir = path.dirname(params.outfile)
	pTsvMerge.output = 'outfile:file:{}'.format(path.basename(params.outfile))
else:
	pPAdjust           = Proc(desc = 'Adjust p values')
	pPAdjust.depends   = pTsvMerge
	pPAdjust.input     = 'infile:file'
	pPAdjust.exdir     = path.dirname(params.outfile)
	pPAdjust.args.pval = params.pval
	pPAdjust.args.fdr  = 'fdr' if params.padj is True else params.padj
	pPAdjust.output    = 'outfile:file:{}'.format(path.basename(params.outfile))
	pPAdjust.lang      = params.Rscript
	pPAdjust.script    = 'file:scripts/ceQTL-pAdjust.R'

PyPPL({
	'default': {'ppldir': params.ppldir, 'forks': params.forks},
	'_log': {'lvldiff': '-DEBUG'}
}).start(pExpmat, pMsmat).run()



