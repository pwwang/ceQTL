"""AtSnp analysis"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect, pTsv
from bioprocs.tfbs import pAtSnp

def pipeline(opts):
    """Construct the pipeline"""

    pTsv.input = [opts.infile]
    pTsv.args.row = 'lambda row: float(row[%r]) < %f' % (opts.pcol, opts.cut)

    pTsvColSelect.depends = pTsv
    pTsvColSelect.args.inopts.cnames = True
    pTsvColSelect.args.cols = [opts.snpcol, opts.tfcol]

    pAtSnp.depends = pTsvColSelect
    pAtSnp.input = lambda ch: ch.cbind(opts.snpbed)
    pAtSnp.args.tflist = opts.tflist
    pAtSnp.args.tfmotifs = opts.motifdb
    pAtSnp.args.plot = False
    pAtSnp.args.nthread = opts.nthread
    pAtSnp.output = ['outfile:file:%s' % Path(opts.outfile).name,
                     'outdir:dir:{{i.infile | fn2}}.atsnp']
    pAtSnp.config.export_dir = Path(opts.outfile).parent
    pAtSnp.config.export_part = 'outfile'

    return pTsv

def main(opts):
    """Main function"""

    PyPPL().start(pipeline(opts)).run()
