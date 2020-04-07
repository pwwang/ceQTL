"""Implementation of ceQTL-tools rs2bed"""
# pylint: disable=assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.snp import pRs2Bed

def pipeline(opts):
    """Construct the pipeline"""
    pRs2Bed.input = [opts.rsfile]
    pRs2Bed.output = 'outfile:file:%s' % Path(opts.outfile).name
    pRs2Bed.args.tool = 'local'
    pRs2Bed.args.inopts = opts.inopts
    pRs2Bed.args.snpcol = opts.snpcol
    pRs2Bed.args.sortby = False
    pRs2Bed.args.ncol = 8
    pRs2Bed.config.export_dir = Path(opts.outfile).parent
    return pRs2Bed

def main(opts):
    """Main function"""
    PyPPL().start(pipeline(opts)).run()
