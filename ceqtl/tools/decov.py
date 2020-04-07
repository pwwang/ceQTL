"""Adjust covariates from expression"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect, pTsvCbind, pMatrixR
from bioprocs.stats import pDeCov

def pipeline(opts):
    """Construct the pipeline"""

    pTsvColSelect.input = opts.covfiles
    pTsvColSelect.args.cols = opts.vars

    pTsvCbind.depends = pTsvColSelect
    pTsvCbind.input = lambda ch: [ch.flatten()]
    pTsvCbind.args.fn2cname = 'function(fn, cnames) cnames'

    pDeCov.depends = pTsvCbind
    pDeCov.tag = opts.tool
    pDeCov.input = lambda ch: ch.insert(0, opts.expr)
    pDeCov.args.nthread = opts.nthread
    pDeCov.args.tool = opts.tool
    pDeCov.args.nk = 0
    pDeCov.lang = opts.peerR

    # mean-centered
    pMatrixR.depends = pDeCov
    pMatrixR.output = 'outfile:file:%s' % Path(opts.outfile).name
    pMatrixR.config.export_dir = Path(opts.outfile).parent
    pMatrixR.args.code = """
    mat = t(scale(t(mat), scale = FALSE))
    """

    return pTsvColSelect

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
