"""Plot ROC curves for different measurements"""

# pylint: disable=invalid-name,assigning-non-slot
import math
from pathlib import Path
from pyppl import PyPPL
from bioprocs.utils import shell2 as shell
from bioprocs.stats import pMediation, pAdjust, pModeration
from bioprocs.tsv import (pTsvHeader,
                          pTsvJoin,
                          pTsvColSelect,
                          pTranspose,
                          pTsvCbind,
                          pTsvSplit,
                          pTsv)
from bioprocs.common import pSort
from procs import pTFTSG2MedCases

def common_samples(gtype, expr, transpose=False):
    """Select common samples between genotype and expression data"""
    pTsvHeaderGT = pTsvHeader.copy()
    pTsvHeaderGT.input = [gtype]

    pTsvHeaderExpr = pTsvHeader.copy()
    pTsvHeaderExpr.input = [expr]

    pSortForJoin = pSort.copy()
    pSortForJoin.depends = pTsvHeaderGT, pTsvHeaderExpr
    pSortForJoin.input = lambda ch1, ch2: ch1.rbind(ch2)

    pTsvJoinHeader = pTsvJoin.copy()
    pTsvJoinHeader.depends = pSortForJoin
    pTsvJoinHeader.input = lambda ch: [ch.flatten()]
    pTsvJoinHeader.args.inopts.cnames = False
    pTsvJoinHeader.args.helper = [
        'rnames_written = False',
        'def write(writer, r1, r2):',
        '   global rnames_written',
        '   if not rnames_written:',
        '       writer.write(["ID"])',
        '       writer.write(["ROWNAME"])',
        '       rnames_written = True',
        '   writer.write(r1)',
    ]
    pTsvJoinHeader.args.do = 'lambda writer, r1, r2: write(writer, r1, r2)'

    pTsvColSelectCommonSamples = pTsvColSelect.copy()
    pTsvColSelectCommonSamples.input = 'infile:file, colfile:file'
    pTsvColSelectCommonSamples.depends = pTsvJoinHeader
    pTsvColSelectCommonSamples.input = lambda ch: ch.rep_row(2).insert(0, [gtype, expr])

    if transpose:
        pTransposeInput = pTranspose.copy()
        pTransposeInput.depends = pTsvColSelectCommonSamples
        return [pTsvHeaderGT, pTsvHeaderExpr], [pTransposeInput]

    return [pTsvHeaderGT, pTsvHeaderExpr], [pTsvColSelectCommonSamples]

def nosplit(opts):
    """Pipeline for no split job"""

    med = 'med' in opts.type.lower()
    starts, ends = common_samples(opts.gtype, opts.expr, True)

    pTFTSG2MedCases.input = opts.tft, opts.snpgene
    starts.append(pTFTSG2MedCases)

    pTsvCbind.depends = ends
    pTsvCbind.input = lambda *chs: [sum((ch.flatten() for ch in chs), [])]
    pTsvCbind.args.fill = False
    pTsvCbind.args.fn2cname = 'function(fn, cnames) cnames'

    pMed = pMediation if med else pModeration
    pMed.depends = pTsvCbind, pTFTSG2MedCases
    pMed.args.pval = opts.pcut
    pMed.args.plot = False
    pMed.args.nthread = opts.nthread
    pMed.output = ('outfile:file:%s, '
                   'outdir:dir:{{i.infile | fn}}.%s') % (
                       Path(opts.outfile).name,
                       'mediation' if med else 'moderation'
                   )
    pMed.config.export_dir = Path(opts.outfile).parent
    pMed.config.export_part = 'outfile'
    return starts

def splitjob(opts):
    """Pipeline for split job"""

    med = 'med' in opts.type.lower()
    starts, ends = common_samples(opts.gtype, opts.expr)

    n_sg_pair = shell.wc_l(opts.snpgene).split()[0]
    # sort the snpgene file by gene for splitting by gene later on
    pSortSG = pSort.copy()
    pSortSG.input = [opts.snpgene]
    pSortSG.args.inopts.skip = 0
    pSortSG.args.params.k = 2
    starts.append(pSortSG)

    pTsvSplit.depends = pSortSG
    pTsvSplit.args.inopts.cnames = False
    pTsvSplit.args.by = math.ceil(float(n_sg_pair)/float(opts.njobs))

    pSortSGBySNP = pSortSG.copy()
    pSortSGBySNP.depends = pTsvSplit
    pSortSGBySNP.input = lambda ch: ch.expand()
    pSortSGBySNP.args.params.k = 1

    # sort genotype file to split
    pSortGT = pSort.copy()
    pSortGT.depends = ends
    pSortGT.input = lambda ch: ch.row_at(0)
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

    pTransposeExpr = pTranspose.copy()
    pTransposeExpr.depends = ends
    pTransposeExpr.input = lambda ch: ch.row_at(1)
    starts.append(pTransposeExpr)

    # sort tft to split
    pSortTFT = pSortSG.copy()
    pSortTFT.input = [opts.tft]
    starts.append(pSortTFT)

    pTFTSplit = pTsvJoin.copy()
    pTFTSplit.depends = pSortTFT, pTsvSplit
    pTFTSplit.input = lambda ch1, ch2: ch2.expand().insert(0, ch1).map(list)
    pTFTSplit.args.inopts.cnames = False
    pTFTSplit.args.match = 'lambda r1, r2: compare(r1[1], r2[1])'
    pTFTSplit.args.do = 'lambda writer, r1, r2: writer.write(r1)'

    pTFTSG2MedCases.depends = pTFTSplit, pTsvSplit
    pTFTSG2MedCases.input = lambda ch1, ch2: ch2.expand().insert(0, ch1)

    pTsvCbind.depends = pTranspose, pTransposeExpr
    pTsvCbind.input = lambda ch1, ch2: ch1.cbind(ch2).map(list)
    pTsvCbind.args.fill = False
    pTsvCbind.args.fn2cname = 'function(fn, cnames) cnames'

    pMed = pMediation if med else pModeration
    pMed.depends = pTsvCbind, pTFTSG2MedCases
    pMed.args.plot = False
    pMed.args.fdr = False
    pMed.args.pval = 1.1

    pAdjust.depends = pMed
    pAdjust.input = lambda ch: [ch.outfile.flatten()]
    pAdjust.args.method = 'BH'
    pAdjust.args.pcol = 'Pval'

    # apply pcut
    pTsv.depends = pAdjust
    pTsv.args.inopts.cnames = True
    pTsv.args.row = 'lambda row: float(row.Pval) < %f %s' % (opts.pcut,
                                                           'and float(row.PropMed) > 0'
                                                           if med else '')

    pTsv.output = 'outfile:file:%s' % Path(opts.outfile).name
    pTsv.config.export_dir = Path(opts.outfile).parent

    return starts

def pipeline(opts):
    """Construct the pipeline"""
    if opts.njobs == 1:
        return nosplit(opts)
    return splitjob(opts)

def main(opts):
    """Main function"""

    PyPPL(forks=opts.njobs).start(pipeline(opts)).run(opts.runner)
