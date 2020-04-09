"""Pathway enrichment analysis"""

# pylint: disable=invalid-name,assigning-non-slot
from pyppl import PyPPL
from bioprocs.tsv import pTsv
from bioprocs.gsea import pEnrichr

def pipeline(opts):
    """Construct the pipeline"""

    pTsv.input = [opts.infile]
    pTsv.args.row = 'lambda row: float(row[%r]) < %f' % (opts.pcol, opts.cut)

    pEnrichr.depends = pTsv
    pEnrichr.args.genecol = int(opts.gcol) if opts.gcol.isdigit() else opts.gcol
    pEnrichr.args.libs = ','.join(opts.dbs)
    pEnrichr.config.export_dir = opts.outdir

    return pTsv

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
