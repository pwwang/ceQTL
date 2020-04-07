"""Implementation of ceQTL-tools simulate"""
# pylint: disable=assigning-non-slot,invalid-name

from pyppl import PyPPL, Proc
from bioprocs.rnaseq import pExprSimulate
from bioprocs.plink import pPlinkSimulate, pPlink2GTMat
from bioprocs.tsv import (pTsvReplaceHeader, pTsv, pTsvSample, pMatrixR,
                          pTsvColSelect)
from bioprocs.stats import pDiffCorr
from bioprocs.bed import pBedIntersect, pBedFlank
from ..procs import pTFT2GeneGroups, pSGPair2Beds

def pipeline(opts):
    """Construct the pipeline"""
    ## Simulate SNPs ange genotypes
    pPlinkSimulate.input = [opts.seed]
    pPlinkSimulate.args.ncases = int(opts.nsamples/2)
    pPlinkSimulate.args.nctrls = opts.nsamples - int(opts.nsamples/2)
    pPlinkSimulate.args.nsnps = opts.nsnps
    pPlinkSimulate.args.minfreq = 0.05

    pPlink2GTMat.depends = pPlinkSimulate
    pPlink2GTMat.args.samid = 'iid'
    pPlink2GTMat.args.snpid = 'raw'

    pTsvReplaceHeader.depends = pPlink2GTMat
    pTsvReplaceHeader.args.cnames = ['ID'] + ['Sample%s' % (i+1)
                                              for i in range(opts.nsamples)]
    pTsvReplaceHeader.output = 'outfile:file:gtype_%s.txt' % opts.seed
    pTsvReplaceHeader.input = lambda ch: ch.col_at(0)
    pTsvReplaceHeader.config.export_dir = opts.outdir

    ## Simulate expression values
    pExprSimulate.input = [opts.seed]
    pExprSimulate.args.nsamples = opts.nsamples
    pExprSimulate.args.ngenes = opts.ngenes + opts.ntfs

    # Change last opts.ntfs to TFx
    pTsv.depends = pExprSimulate
    # TF ID startswith 0 as SNP ID does
    pTsv.args.row = ('lambda row: row '
                     'if int(row[0][4:]) <= {0} '
                     'else ["TF" + str(int(row[0][4:]) - {0} - 1)] + '
                     'row[1:]'.format(opts.ngenes))
    pTsv.output = 'outfile:file:expr_%s.txt' % opts.seed
    pTsv.config.export_dir = opts.outdir

    ## Simulate TF/Snp-gene pairs
    pAllPairs = Proc(desc='Generate files with all tf/snp-gene pairs')
    pAllPairs.input = {
        'ngenes,nother,pergene,name': [
            (opts.ngenes, opts.ntfs, opts.tfpergene, 'TF'),
            (opts.ngenes, opts.nsnps, opts.snppergene, 'SNP_'),
        ]
    }
    pAllPairs.output = ('outfile:file:{{i.name | .rstrip: "_" '
                        '                      | @append: "-Gene.txt"}}')
    pAllPairs.forks = opts.nthread
    pAllPairs.script = '''
    for iother in $(seq 0 {{i.nother | @minus: 1}}); do
        for igene in $(seq 1 {{i.ngenes}}); do
            echo -e "{{i.name}}$iother\\tGene$igene" >> {{o.outfile | quote}}
        done
    done
    '''

    # Turn genotype matrix into binary for each SNP
    pMatrixR.depends = pTsvReplaceHeader
    pMatrixR.args.code = '''
    mat[mat>1] = 1
    mat = t(mat)
    '''

    pTFTSample = pTsvSample.copy()
    pTFTSample.depends = pAllPairs
    pTFTSample.input = lambda ch: ch.row_at(0)
    pTFTSample.output = 'outfile:file:tft_{{args.seed}}.txt'
    pTFTSample.args.seed = opts.seed
    pTFTSample.args.inopts.cnames = False
    pTFTSample.args.n = opts.ngenes * opts.tfpergene
    pTFTSample.config.export_dir = opts.outdir

    pSGPair2Beds.depends = pAllPairs
    pSGPair2Beds.input = lambda ch: ch.row_at(1)
    pSGPair2Beds.config.export_dir = opts.outdir
    pSGPair2Beds.args.snppergene = opts.snppergene
    pSGPair2Beds.args.dist = 1000

    pBedFlank.depends = pSGPair2Beds
    pBedFlank.input = lambda ch: ch.genefile
    pBedFlank.args.extend = True
    pBedFlank.args.params.b = 1000

    pBedIntersect.depends = pSGPair2Beds, pBedFlank
    pBedIntersect.input = lambda ch1, ch2: ch1.snpfile.cbind(ch2)
    pBedIntersect.args.params.wa = True
    pBedIntersect.args.params.wb = True

    pTsvColSelect.depends = pBedIntersect
    pTsvColSelect.args.inopts.cnames = False
    pTsvColSelect.args.cols = [3, 9]
    pTsvColSelect.output = 'outfile:file:snpgene_%s.txt' % opts.seed
    pTsvColSelect.config.export_dir = opts.outdir

    pTFT2GeneGroups.depends = pTFTSample

    pDiffCorr.depends = pTsv, pMatrixR, pTsvColSelect, pTFT2GeneGroups
    pDiffCorr.config.export_dir = opts.outdir
    pDiffCorr.args.nthread = opts.nthread

    return [pPlinkSimulate,
            pExprSimulate,
            pAllPairs]

def main(opts):
    """Main function"""
    PyPPL().start(pipeline(opts)).run()
