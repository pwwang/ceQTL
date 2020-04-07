"""Plot ROC curves for different measurements"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
import cmdy
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect, pTsvHeader, pTsv, pTsvCbind, pTsvReplaceHeader
from bioprocs.common import pSort, pShell
from bioprocs.stats import pHypergeom

def colselect(infile, col, cut):
    """Select the columns from input file"""
    pTsvColSelect1 = pTsvColSelect.copy()
    if str(col).isdigit():
        pTsvColSelect1.input = [infile]
        pTsvColSelect1.args.cols = [0, int(col)]
        start = end = pTsvColSelect1
    else:
        pTsvHeader1 = pTsvHeader.copy()
        pTsvHeader1.input = [infile]
        start = pTsvHeader1

        pShell1 = pShell.copy()
        pShell1.depends = pTsvHeader1
        pShell1.args.cmd = """
        head -1 {{i.infile | quote}} > {{o.outfile | quote}}
        echo "%s" >> {{o.outfile | quote}}
        """ % col

        pTsvColSelect1.depends = pShell1
        pTsvColSelect1.input = 'infile:file, colfile:file'
        pTsvColSelect1.input = lambda ch: ch.insert(0, infile)
        start = pTsvHeader1
        end = pTsvColSelect1

    if cut and cut < 1:
        pTsv1 = pTsv.copy()
        pTsv1.depends = pTsvColSelect1
        pTsv1.args.row = 'lambda row: float(row[%r]) < %f' % (col, cut)
        end = pTsv1
    elif cut and cut > 1:
        pSort1 = pSort.copy()
        pSort1.depends = pTsvColSelect1
        pSort1.args.inopts.skip = 1
        pSort1.args.params.k = '2g'

        pShell2 = pShell.copy()
        pShell2.depends = pSort1
        pShell2.args.cmd = "head -n %d {{i.infile | quote}} > {{o.outfile | quote}}" % (int(cut)+1)
        end = pShell2

    return start, end

def pipeline(opts):
    """Construct the pipeline"""

    start, end = colselect(opts.infile, opts.col, opts.cut)

    # duplicate gold standard columns
    pTsvGold = pTsv.copy()
    pTsvGold.input = [opts.gold]
    pTsvGold.args.inopts.cnames = False
    pTsvGold.args.row = 'lambda row: [row[0], 1]'
    start = [pTsvGold, start]

    # add header
    pTsvReplaceHeader.depends = pTsvGold
    pTsvReplaceHeader.args.inopts.cnames = False
    pTsvReplaceHeader.args.cnames = ['ROWNAME', 'GOLD']

    pTsvCbind.depends = pTsvReplaceHeader, end
    pTsvCbind.input = (lambda *chs: [[ch.get() for ch in chs]])

    # prepare file for hypergeometric test
    # replace NA with 0
    # replace pvalues with 1 (presence)
    pPrepHG = pTsv.copy()
    pPrepHG.depends = pTsvCbind
    pPrepHG.args.row = ('lambda row: row.__setitem__(1, int(row[1] == "1")) '
                        'or row.__setitem__(2, int(row[2] != "NA"))')

    pHypergeom.depends = pPrepHG
    pHypergeom.args.intype = 'raw'
    if str(opts.bign).isdigit():
        pHypergeom.args.N = int(opts.bign)
    else:
        bign = cmdy.wc(l = opts.bign).strip().split()[0]
        pHypergeom.args.N = int(bign)

    return start

def main(opts):
    """Main function"""

    PyPPL(forks=opts.nthread).start(pipeline(opts)).run()
