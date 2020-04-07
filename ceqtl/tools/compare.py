"""Rank compare of two measurements (i.e. ceQTL, eQTL, median, etc.)"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect, pTsvJoin, pTsvHeader, pTsv
from bioprocs.plot import pScatterCompare
from bioprocs.common import pSort, pShell

def colselect(infile, col, cut):
    """Select the columns from input file"""
    tag = Path(infile).stem.split(".")[0]
    pTsvColSelect1 = pTsvColSelect.copy(tag=tag)
    if str(col).isdigit():
        pTsvColSelect1.input = [infile]
        pTsvColSelect1.args.cols = [0, int(col)]
        start = pTsvColSelect1
    else:
        pTsvHeader1 = pTsvHeader.copy(tag=tag)
        pTsvHeader1.input = [infile]
        start = pTsvHeader1

        pShell1 = pShell.copy(tag=tag)
        pShell1.depends = pTsvHeader1
        pShell1.args.cmd = """
        head -1 {{i.infile | quote}} > {{o.outfile | quote}}
        echo "%s" >> {{o.outfile | quote}}
        """ % col

        pTsvColSelect1.depends = pShell1
        pTsvColSelect1.input = lambda ch: ch.insert(0, infile)

    pTsv1 = pTsv.copy(tag=tag)
    pTsv1.depends = pTsvColSelect1
    pTsv1.args.row = "lambda row: float(row[1]) < %f" % cut
    return [start], [pTsv1]

def pipeline(opts):
    """Construct the pipeline"""

    starts1, ends1 = colselect(opts.in1, opts.col1, opts.cut1)
    starts2, ends2 = colselect(opts.in2, opts.col2, opts.cut2)
    starts = starts1 + starts2
    ends = ends1 + ends2

    pSort.depends = ends
    pSort.input = lambda ch1, ch2: ch1.flatten() + ch2.flatten()
    pSort.args.inopts.skip = 1

    pTsvJoin.depends = pSort
    pTsvJoin.input = lambda ch: [ch.flatten()]
    pTsvJoin.args.outopts.cnames = ["ROWNAME", Path(opts.in1).stem, Path(opts.in2).stem]
    pTsvJoin.args.inopts.cnames = True
    pTsvJoin.args.do = "lambda writer, r1, r2: writer.write(r1[:2] + [r2[1]])"

    pScatterCompare.depends = pTsvJoin
    pScatterCompare.args.stacked = False
    pScatterCompare.args.tsform = 'function(x) order(-x)'

    return starts

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
