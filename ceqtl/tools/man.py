"""Manhattan plot"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect, pTsvJoin, pTsvAggregate
from bioprocs.plot import pManhattan2
from bioprocs.common import pSort

def pipeline(opts):
    """Construct the pipeline"""

    pTsvColSelect.input = [opts.infile]
    pTsvColSelect.args.cols = ["Case", opts.col]

    pTsvAggregate.depends = pTsvColSelect
    pTsvAggregate.args.sort = False
    pTsvAggregate.args.aggrs.Padj = "$min:1"
    pTsvAggregate.args.on = "Case"

    pSortInfile = pSort.copy()
    pSortInfile.depends = pTsvAggregate
    pSortInfile.args.inopts.skip = 1

    pSortPos = pSort.copy()
    pSortPos.input = [opts.pos]
    pSortPos.args.params.k = 4

    pTsvJoin.depends = pSortInfile, pSortPos
    pTsvJoin.input = lambda ch1, ch2: [ch1.flatten() + ch2.flatten()]
    pTsvJoin.args.inopts.cnames = True, False
    pTsvJoin.args.outopts.cnames = ['SNP', 'Chr', 'Position', 'Padj']
    pTsvJoin.args.match = 'lambda r1, r2: compare(r1[0], r2[3])'
    pTsvJoin.args.do = ('lambda writer, r1, r2: writer.write('
                        '   [r2[3], r2[0][3:] if r2[0][:3] == "chr" else r2[0], r2[2], r1[1]]'
                        ')')

    if not opts.params:
        opts.params = {}
    opts.params.setdefault("ylab", "-log10(Padj)")
    pManhattan2.depends = pTsvJoin
    pManhattan2.args.hifile = opts.hilight
    pManhattan2.args.params = opts.params
    pManhattan2.args.devpars.height = 1200
    pManhattan2.args.devpars.res = 150
    pManhattan2.output = 'outfile:file:%s' % Path(opts.outfile).name
    pManhattan2.config.export_dir = Path(opts.outfile).parent

    return pTsvColSelect, pSortPos

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
