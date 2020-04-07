"""Q-Q plot for ceQTL results"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect
from bioprocs.plot import pQQ

def pipeline(opts):
    """Construct the pipeline"""

    pTsvColSelect.input = [opts.infile]
    pTsvColSelect.args.cols = ['Case', opts.col]

    pQQ.depends = pTsvColSelect
    pQQ.args.inopts.rnames = True
    pQQ.args.tsform = 'function(x) -log10(x)'
    pQQ.args.ggs.ylab = {0: "-log10(%s)" % opts.col}
    pQQ.output = 'outfile:file:%s' % Path(opts.outfile).name
    pQQ.config.export_dir = Path(opts.outfile).parent

    return pTsvColSelect

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
