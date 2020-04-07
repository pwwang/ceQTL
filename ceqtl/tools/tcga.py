"""Implementation of ceQTL-tools tcga"""
# pylint: disable=assigning-non-slot
from os import path
from pyppl import PyPPL
from bioprocs.tcga import pTCGADownload, pSample2SubmitterID, pGtFiles2Mat
from bioprocs.common import pShell
from bioprocs.web import pDownloadGet
from bioprocs.tsv import pTsvColSelect
from bioprocs.gene import pGeneNameNorm

def pipeline(opts):
    """Construct the pipeline"""
    outdir = path.join(opts.outdir, opts.cancer)

    pTCGADownload.input = [opts.snpmani]
    pTCGADownload.args.nthread = opts.nthread
    pTCGADownload.args.token = opts.token
    pTCGADownload.config.export_dir = outdir
    pTCGADownload.cache = 'export'

    pSample2SubmitterID.depends = pTCGADownload
    pSample2SubmitterID.input = lambda ch: ch.cbind(opts.snpmeta)
    pSample2SubmitterID.args.len = 16

    # remove normal-like
    pShell.depends = pSample2SubmitterID
    pShell.args.cmd = '''
    mkdir {{o.outfile}}
    for gtfile in {{i.infile}}/*.txt; do
        if [[ $(basename $gtfile | cut -c14) == "0" ]]; then
            ln -s $gtfile {{o.outfile}}/
        fi
    done
    '''

    pGtFiles2Mat.depends = pShell
    pGtFiles2Mat.input = lambda ch: [ch.expand(pattern='*.txt').flatten()]
    pGtFiles2Mat.config.echo_jobs = [0]
    pGtFiles2Mat.config.echo_types = 'all'
    pGtFiles2Mat.config.export_dir = outdir
    pGtFiles2Mat.output = 'outfile:file:TCGA-%s.gt.txt' % opts.cancer

    pDownloadGet.input = [
        'https://gdc.xenahubs.net/download/'
        'TCGA-%s.htseq_counts.tsv.gz' % opts.cancer,
        'https://gdc.xenahubs.net/download/'
        'TCGA-%s.GDC_phenotype.tsv.gz' % opts.cancer,
        'https://gdc.xenahubs.net/download/'
        'TCGA-%s.survival.tsv.gz' % opts.cancer
    ]
    pDownloadGet.config.export_dir = outdir
    pDownloadGet.config.export_part = ['*phenotype.tsv.gz',
                                       '*.survival.tsv.gz']

    # remove normal-like
    pTsvColSelect.depends = pDownloadGet
    pTsvColSelect.input = lambda ch: ch.row_at(0)
    pTsvColSelect.args.cols = ('lambda cnames: [name for name in cnames '
                               'if name[-3:-2] == "0" or name[:5] != "TCGA-"]')

    # convert ENSG to gene symbols
    pGeneNameNorm.depends = pTsvColSelect
    pGeneNameNorm.output = 'outfile:file:TCGA-%s.expr.txt' % opts.cancer
    pGeneNameNorm.args.inopts.cnames = True
    pGeneNameNorm.args.frm = 'ensembl.gene'
    pGeneNameNorm.args.notfound = 'skip'
    pGeneNameNorm.config.export_dir = outdir

    return pTCGADownload, pDownloadGet

def main(opts):
    """Main function"""
    PyPPL().start(pipeline(opts)).run()
