"""Implementation of ceQTL-tools genotype-qc"""
# pylint: disable=invalid-name,assigning-non-slot

from pathlib import Path
from pyppl import PyPPL
from bioprocs.plink import (pGTMat2Plink,
                            pPlinkMiss,
                            pPlinkFreq,
                            pPlinkHWE,
                            pPlinkIBD,
                            pPlinkRemove,
                            pPlink2GTMat,
                            pPlinkStats)

def plink_from_gtmat(opts):
    """Convert genotype matrix to plink format"""
    pGTMat2Plink.input = [opts.gtmat]
    pGTMat2Plink.args.snpbed = opts.snpbed

    return pGTMat2Plink

def plink_to_gtmat(opts, indir):
    """Convert plink back to GT matrix"""
    pPlinkStats.input = [indir]
    pPlinkStats.args.params['het'] = False
    pPlinkStats.args.params['check-sex'] = False
    pPlinkStats.args.cutoff['hardy.hwe'] = opts.hwe
    pPlinkStats.args.cutoff['missing.sample'] = opts.samplecr
    pPlinkStats.args.cutoff['missing.snp'] = opts.snpcr
    pPlinkStats.args.cutoff['freq'] = opts.maf
    pPlinkStats.args.plot['het'] = False
    pPlinkStats.args.plot['hardy.mingt'] = False
    pPlinkStats.config.export_dir = opts.outdir

    pPlink2GTMat.input = [indir]
    pPlink2GTMat.output = 'outfile:file:%s, outsnp:file:%s.snp.bed' % (
        Path(opts.gtmat).name, Path(opts.gtmat).stem.split('.')[0])
    pPlink2GTMat.args.samid = 'iid'
    pPlink2GTMat.args.snpid = '{rs}'
    pPlink2GTMat.args.refallele = opts.snpbed
    pPlink2GTMat.config.export_dir = opts.outdir
    return pPlinkStats, pPlink2GTMat

def pipeline(opts, tag, indir):
    """Construct the pipeline"""

    # processes are not supposed to be reused
    # so we copy each process

    # SNP & sample call rate
    pPlinkMissQC = pPlinkMiss.copy(tag=tag)
    pPlinkMissQC.input = [indir]
    pPlinkMissQC.args.samplecr = opts.samplecr
    pPlinkMissQC.args.snpcr = opts.snpcr
    pPlinkMissQC.args.plot = False

    pPlinkRemoveMiss = pPlinkRemove.copy(tag=tag)
    pPlinkRemoveMiss.depends = pPlinkMissQC
    pPlinkRemoveMiss.input = lambda ch: ch.insert(0, indir)

    pPlinkFreqQC = pPlinkFreq.copy(tag=tag)
    pPlinkFreqQC.depends = pPlinkRemoveMiss
    pPlinkFreqQC.args.cutoff = opts.maf
    pPlinkFreqQC.args.plot = False

    pPlinkRemoveFreq = pPlinkRemove.copy(tag=tag)
    pPlinkRemoveFreq.depends = pPlinkRemoveMiss, pPlinkFreqQC

    pPlinkHWEQC = pPlinkHWE.copy(tag=tag)
    pPlinkHWEQC.depends = pPlinkRemoveFreq
    pPlinkHWEQC.args.cutoff = opts.hwe
    pPlinkHWEQC.args.plot = False

    pPlinkRemoveHWE = pPlinkRemove.copy(tag=tag)
    pPlinkRemoveHWE.depends = pPlinkRemoveFreq, pPlinkHWEQC

    pPlinkIBDQC = pPlinkIBD.copy(tag=tag)
    pPlinkIBDQC.depends = pPlinkRemoveHWE
    pPlinkIBDQC.args.plot = False
    pPlinkIBDQC.args.seed = 8525

    pPlinkRemoveIBD = pPlinkRemove.copy(tag=tag)
    pPlinkRemoveIBD.depends = pPlinkRemoveHWE, pPlinkIBDQC
    pPlinkRemoveIBD.output = 'outdir:dir:%s.plink' % Path(opts.gtmat).stem

    return pPlinkMissQC, pPlinkRemoveIBD

def main(opts):
    """Main function"""
    PyPPL(name='Convert genotype matrix to plink format') \
        .start(plink_from_gtmat(opts)) \
        .run()

    indir = pGTMat2Plink.channel.get(0)

    for i in range(opts.iter):
        start, end = pipeline(opts, tag='iter%s' % (i + 1), indir=indir)
        if i == opts.iter - 1:
            end.config.export_dir = opts.outdir

        PyPPL(name="Genotype data quality control, iteration %s" % (i + 1)) \
            .start(start) \
            .run()

        indir = end.channel.get(0)

    PyPPL(name='Convert plink back to genotype matrix') \
        .start(plink_to_gtmat(opts, indir)) \
        .run()
