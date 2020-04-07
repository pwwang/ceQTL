"""Whole genome scan using given motifs and mutations"""
# pylint: disable=assigning-non-slot,redefined-outer-name,invalid-name

from pyppl import PyPPL

from bioprocs.seq import pPromoters, pSeqMutate
from bioprocs.bed import pBedGetfasta
from bioprocs.tfbs import pMotifScan
from bioprocs.common import pShell

def pipeline(opts):
    """Main function"""

    pPromoters.input = [opts.gfile]
    pPromoters.args.region.up = opts.region
    pPromoters.args.region.withbody = True
    pPromoters.args.genecol = opts.gcol
    pPromoters.args.inopts = opts.inopts

    pBedGetfasta.depends = pPromoters
    pBedGetfasta.args.params['name+'] = True

    pMotifScan.input = lambda ch: ch.insert(0, opts.tflist)
    pMotifScan.args.pval = opts.pcut
    pMotifScan.args.nthread = opts.nthread
    pMotifScan.config.export_dir = opts.outdir

    if opts.mut:
        pSeqMutate.depends = pBedGetfasta
        pSeqMutate.dirsig = False
        pSeqMutate.input = lambda ch: ch.cbind(opts.mut)
        pSeqMutate.args.nthread = opts.nthread
        pMotifScan.depends = pSeqMutate
    else:
        pMotifScan.depends = pBedGetfasta

    pShell.depends = pMotifScan
    pShell.input = lambda ch: ch.outfile
    pShell.output = 'outfile:file:{{i.infile | stem2}}.tf-gene.txt'
    pShell.args.cmd = ('grep -v "^#" {{i.infile | quote}} '
                       '| cut -f4 '
                       '| sed "s/::/\\t/" > {{o.outfile | quote}}')
    pShell.config.export_dir = opts.outdir

    return pPromoters

def main(opts):
    """Main function"""
    PyPPL().start(pipeline(opts)).run()
