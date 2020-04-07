"""Manhattan plot"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
import cmdy
from pyppl import PyPPL
from bioprocs.tsv import pTsvJoin, pTsvAggregate
from bioprocs.common import pSort

def pipeline(opts):
    """Construct the pipeline"""

    # get the column index of the column name
    colidx = cmdy.head(opts.infile, n=1, _pipe=True) | \
        cmdy.tr("\t", "\n", _pipe=True) | \
        cmdy.nl(_pipe=True) | \
        cmdy.grep(w=opts.on)
    colidx = colidx.strip().split()[0]

    # sort by that column, make sure the first element we remain
    # is the corresponding record
    pSort.input = [opts.infile]
    pSort.args.inopts.skip = 1
    pSort.args.params.k = '%sg' % colidx

    pTsvAggregate.depends = pSort
    pTsvAggregate.args.sort = False
    pTsvAggregate.args.aggrs.Padj_aggr = "$%s:%s" % (opts.method, opts.on)
    #pTsvAggregate.args.on = "Case"
    pTsvAggregate.args.origin = "keep0"
    pTsvAggregate.output = "outfile:file:%s" % Path(opts.outfile).name
    pTsvAggregate.config.export_dir = Path(opts.outfile).parent

    return pSort

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
