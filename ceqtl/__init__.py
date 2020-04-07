"""ceQTL pipeline"""
# pylint: disable=assigning-non-slot
import math
from pathlib import Path
from diot import Diot
from pyppl import PyPPL
from bioprocs import params

__version__ = "0.0.1"

params._desc = __doc__
params.expr.required = True
params.expr.desc = ('The expression matrix, '
                    'with genes as rows and samples as columns.')
params.gtype.required = True
params.gtype.desc = ('The genotype matrix, '
                     'with SNPs as rows and samples as columns.')
params.tft.required = True
params.tft.desc = ('The TF target matrix, '
                   'with target genes as rows and TFs as columns.')

params.snpgene.required = True
params.snpgene.desc = [
    'The SNP-gene pairs to limit the number of regressions.',
    'This is usually decided by the distance of the SNP to the gene.',
    'You can use `ceqtl-tools snpgene` to generate this file.'
]

params.pcut = 0.05
params.pcut.desc = 'The pvalue cutoff for the trios.'
params.pval = params.pcut

params.njobs = 1
params.njobs.desc = ('Split the cases into different jobs. '
                     'You can distribute the jobs into clusters, for example.')

params.ncores = 1
params.ncores.desc = 'Number of cores to use each job.'
params.nthread = params.ncores

params.runner = 'local'
params.runner.desc = ('The runner for the pipeline. '
                      'See https://pyppl.readthedocs.io/en/latest/runners/')

params.outfile.required = True
params.outfile.desc = "The output file"

def main():
    """Main function"""
    opts = params._parse(dict_wrapper=Diot)
    from bioprocs.stats import pChow, pAdjust
    from bioprocs.tsv import pTranspose, pTsvSplit, pTsvJoin, pTsv
    from bioprocs.common import pSort
    from bioprocs.utils import shell2 as shell
    from procs import pTFT2GeneGroups

    if (opts.njobs == 1):
        pTranspose.input = [opts.gtype]
        pTFT2GeneGroups.input = [opts.tft]

        pChow.depends = pTFT2GeneGroups, pTranspose
        pChow.input = lambda ch1, ch2: ch1.insert(0, opts.expr, ch2.get(),
                                                  opts.snpgene)
        pChow.args.nthread = opts.ncores
        pChow.args.pval = opts.pcut
        pChow.args.plot = False
        pChow.output = ('outfile:file:%s, '
                        'outdir:dir:{{i.infile | fn}}.chow') % Path(opts.outfile).name
        pChow.config.export_dir = Path(opts.outfile).parent

        start_processes = pTranspose, pTFT2GeneGroups
    else:
        n_sg_pair = shell.wc_l(opts.snpgene).split()[0]
        # sort the snpgene file by gene for splitting by gene later on
        pSortSG = pSort.copy()
        pSortSG.input = [opts.snpgene]
        pSortSG.args.inopts.skip = 0
        pSortSG.args.params.k = 2

        pTsvSplit.depends = pSortSG
        pTsvSplit.args.inopts.cnames = False
        pTsvSplit.args.by = math.ceil(float(n_sg_pair)/float(opts.njobs))

        pSortSGBySNP = pSortSG.copy()
        pSortSGBySNP.input = lambda ch: ch.expand()
        pSortSGBySNP.depends = pTsvSplit
        pSortSGBySNP.args.params.k = 1

        # sort genotype file to split
        pSortGT = pSort.copy()
        pSortGT.input = [opts.gtype]
        pSortGT.args.inopts.skip = 1
        pSortGT.args.params.k = 1

        # select gtype type for each set of genes
        pGTSplit = pTsvJoin.copy()
        pGTSplit.depends = pSortGT, pSortSGBySNP
        # size:                  1,  opts.njobs
        pGTSplit.input = lambda ch1, ch2: ch2.insert(0, ch1).map(list)
        pGTSplit.args.inopts.cnames = [True, False]
        pGTSplit.args.outopts.cnames = 0
        pGTSplit.args.match = 'lambda r1, r2: compare(r1[0], r2[0])'
        pGTSplit.args.do = 'lambda writer, r1, r2: writer.write(r1)'

        pTranspose.depends = pGTSplit

        # sort tft to split
        pSortTFT = pSortSG.copy()
        pSortTFT.input = [opts.tft]

        pTFTSplit = pTsvJoin.copy()
        pTFTSplit.depends = pSortTFT, pTsvSplit
        pTFTSplit.input = lambda ch1, ch2: ch2.expand().insert(0, ch1).map(list)
        pTFTSplit.args.inopts.cnames = False
        pTFTSplit.args.match = 'lambda r1, r2: compare(r1[1], r2[1])'
        pTFTSplit.args.do = 'lambda writer, r1, r2: writer.write(r1)'

        pTFT2GeneGroups.depends = pTFTSplit

        pChow.depends = pTranspose, pTsvSplit, pTFT2GeneGroups
        pChow.input = (lambda ch1, ch2, ch3:
                       ch1.cbind(ch2.expand(), ch3).insert(0, opts.expr))
        pChow.args.plot = False
        pChow.args.fdr = False # don't do fdr for single job
        pChow.args.pval = 1.1 # all pvalues to calculate adjusted p

        pAdjust.depends = pChow
        pAdjust.input = lambda ch: [ch.outfile.flatten()]
        pAdjust.args.method = 'BH'
        pAdjust.args.pcol = 'Pval'

        # apply pcut
        pTsv.depends = pAdjust
        pTsv.args.inopts.cnames = True
        pTsv.args.row = 'lambda row: float(row.Pval) < %f' % opts.pcut

        pTsv.output = 'outfile:file:%s' % Path(opts.outfile).name
        pTsv.config.export_dir = Path(opts.outfile).parent

        start_processes = [pSortSG, pSortGT, pSortTFT]

    PyPPL(forks=opts.njobs).start(start_processes).run(opts.runner)

if __name__ == '__main__':
    main()