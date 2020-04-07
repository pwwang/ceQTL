"""Implementation of ceQTL-tools pairsnpgenes"""
# pylint: disable=assigning-non-slot,redefined-outer-name,invalid-name
from pathlib import Path
from pyppl import PyPPL

from bioprocs.seq import pPromoters
from bioprocs.bed import pBedIntersect
from bioprocs.common import pShell

def pipeline(opts):
    """Main function"""

    pPromoters.input = [opts.gfile]
    pPromoters.args.region.up = opts.region
    pPromoters.args.region.withbody = True
    pPromoters.args.genecol = opts.gcol
    pPromoters.args.inopts = opts.inopts

    pBedIntersect.depends = pPromoters
    pBedIntersect.input = lambda ch: ch.cbind(opts.snpbed)
    pBedIntersect.args.wa = True
    pBedIntersect.args.wb = True

    pShell.depends = pBedIntersect
    pShell.output = 'outfile:file:%s' % Path(opts.outfile).name
    pShell.config.export_dir = Path(opts.outfile).parent
    pShell.args.cmd = ("awk 'OFS=\"\\t\"{print $10,$4}' {{i.infile | quote}} "
                       "| sort -k2 -u > {{o.outfile | quote}}")

    return pPromoters

def main(opts):
    """Main function"""
    PyPPL().start(pipeline(opts)).run()
