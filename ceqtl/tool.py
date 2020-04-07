"""A couple of tools to prepare data for ceQTL pipeline"""
from pathlib import Path
from pyparam import commands
from diot import Diot
from bioprocs import params

### Download TCGA data
commands.tcga = [
    'Download TCGA data and prepare them.',
    'Currently only TCGA data supported.',
    'We will download clinic data and expression data from UCSC Xena and '
    'SNP6 array data from TCGA legacy archive to get genotypes.'
    'Note that the expression unit is `log2(count+1)`'
]
commands.tcga.cancer.required = True
commands.tcga.cancer.desc = 'The cancer id from TCGA. For example: LUAD'
commands.tcga.snpmani.required = True
commands.tcga.snpmani.desc = 'The manifest file for SNP genotyping data'
commands.tcga.snpmeta.required = True
commands.tcga.snpmeta.desc = 'The metadata file for SNP genotyping data'
commands.tcga.token.desc = 'Token file for download TCGA files.'

### Genotype data QC
commands['genotype-qc'] = [
    'Quality control for genotype data by:',
    '- [SNP] SNP call rate (called acorss all samples)',
    '- [Sample] Sample call rate (how many SNPs call of all SNPs)',
    '- [SNP] Minor allele frequency',
    '- [SNP] Hardy Weinberg equilibrium test',
    '- [Sample] Identity by state (IBS)'
]
GTQC = commands['genotype-qc']
GTQC.gtmat.required = True
GTQC.gtmat.desc = 'The genotype matrix with rows as SNPs and columns as samples'
GTQC.snpcr = 0.95
GTQC.snpcr.desc = 'SNP call rate cutoff'
GTQC.hwe = 1e-4
GTQC.hwe.desc = 'Hardy Weinberg test cutoff'
GTQC.maf = 0.1
GTQC.maf.desc = 'Minor allele frequency cutoff'
GTQC.samplecr = .9
GTQC.samplecr.desc = 'Sample call rate cutoff'
GTQC.ibs = 0.95
GTQC.ibs.desc = 'Identity by state cutoff'
GTQC.iter = 3
GTQC.iter.desc = 'Number of iterations to do the filtering'
GTQC.snpbed.required = True
GTQC.snpbed.desc = 'The BED6+ file for SNP coordinates'
commands['gtype-qc'] = GTQC

### Expression data QC
commands['expr-qc'] = [
    'Quality control for expression data by:',
    '- Minimum median expression value before transformation',
    '- Minimum median expression value after transformation',
    '- The workflow:',
    '  1. pre-transformation median value filtering',
    '  2. input transformation (i.e convert log2 scale back to normal)',
    '  3. unit conversion (i.e. convert raw counts to tmm)',
    '  4. output transformation (i.e convert to log2 scale)',
    '  5. post-transformation median value filtering',
]
EXPRQC = commands['expr-qc']
EXPRQC.expr.required = True
EXPRQC.expr.desc = ('The expression matrix with rows as genes '
                    'and columns as samples')
EXPRQC.premed = 1e-9
EXPRQC.premed.desc = ('The minimum median expression before transformation.'
                      'Set it to `0` to skip filtering.')
EXPRQC.postmed = 1
EXPRQC.postmed.desc = ('The minimum median expression after transformation.'
                       'Set it to `0` to skip filtering.')
EXPRQC.inform.desc = ('An R function to transform the input data.'
                      'For example, convert log2 scale to normal.')
EXPRQC.inunit = 'count'
EXPRQC.inunit.desc = 'The unit of the input expression values'
EXPRQC.outunit = 'tmm'
EXPRQC.outunit.desc = ('The unit of the output expression values.'
                       'If you don`t want to do unit conversion, '
                       'set it to `False`, then `outtsform` will be ignored.')
EXPRQC.outform = 'log2'
EXPRQC.outform.desc = ('An R function to transform output data. '
                       '`log2` is short for `function(expr) log2(expr+1)`')

### Extract coordinates from SNP ids
commands.rs2bed = 'Convert rs IDs to a bed6+ file'
commands.rs2bed.rsfile.required = True
commands.rs2bed.rsfile.desc = 'The file with rs IDs'
commands.rs2bed.inopts = Diot(delimit="\t", skip=0, comment="#")
commands.rs2bed.inopts.desc = 'Options to read the rsfile'
commands.rs2bed.snpcol = 0
commands.rs2bed.snpcol.desc = ('The column index (0-based) '
                               'where the rs IDs located')
commands.rs2bed.dbsnp = params.dbsnp
commands.rs2bed.dbsnp.show = True
commands.rs2bed.outfile.desc = (
    'The output file. If this is specified with a full path, '
    '`--outdir` will be ignored. Otherwise the outfile will be '
    '`<outdir>/<outfile>`'
)
commands.rs2bed.outfile.callback = (
    lambda opt, ps: opt.set_value(Path(ps.outdir.value).joinpath(opt.value))
    if ps.outdir.value and opt.value and '/' not in opt.value
    else opt.set_value(Path(ps.outdir.value).joinpath(
        Path(ps.rsfile).with_suffix('.snp.bed')
    ))
    if ps.outdir.value and not opt.value
    else None
    if '/' in opt.value
    else 'Either --outdir or --outfile must be specified'
)

### PWM scan for potential TF-gene pairs
commands.pwmscan = 'PWM scan with given motifs and mutations'
commands.pwmscan.tfmotifs = params.tfmotifs
commands.pwmscan.tfmotifs.show = True
commands.pwmscan.tflist = params.tflist
commands.pwmscan.tflist.show = True
commands.pwmscan.gfile.required = True
commands.pwmscan.gfile.desc = 'The file containing genes'
commands.pwmscan.inopts = Diot(delimit="\t", skip=0, comment="#")
commands.pwmscan.inopts.desc = 'Options to read the gfile'
commands.pwmscan.gcol = 0
commands.pwmscan.gcol.desc = ('The column index (0-based) '
                              'where the gene located')
commands.pwmscan.ref = params.ref
commands.pwmscan.ref.show = True
commands.pwmscan.region = 250000
commands.pwmscan.region.desc = 'The flank region of a gene to consider'
commands.pwmscan.mut.desc = ('The mutations in BED6+ format with 7th column '
                             'the reference allele and 8th the alternative')
commands.pwmscan.pcut = 1e-6
commands.pwmscan.pcut.desc = 'The p-value cutoff'

### Generate simulated data
commands.simulate = 'Simulate the data for evaluation.'
commands.simulate.seed = 8525
commands.simulate.seed.desc = 'The seed for the simulation'
commands.simulate.nsamples = 100
commands.simulate.nsamples.desc = 'Number of samples to simulate'
commands.simulate.nsnps = 200
commands.simulate.nsnps.desc = 'Number of snps to simulate'
commands.simulate.ngenes = 1000
commands.simulate.ngenes.desc = 'Number of genes to simulate'
commands.simulate.snppergene = 5
commands.simulate.snppergene.desc = 'Number cis-SNPs per gene on average'
commands.simulate.tfpergene = 2
commands.simulate.tfpergene.desc = 'Number TFs per gene on average'
commands.simulate.ntfs = 20
commands.simulate.ntfs.desc = 'Number of TFs to simulate'

### Add and sort TF-gene pairs
commands.sorttfgenes = ('Add TF-gene pairs to existing ones, '
                        'remove redundant ones, and sort them by genes')
commands.sorttfgenes.origin.required = True
commands.sorttfgenes.origin.desc = ('The original tf-gene pairs '
                                    'with 1st column TF and 2nd column Gene')
commands.sorttfgenes.addition.required = True
commands.sorttfgenes.addition.desc = ('The tf-gene pairs to add, '
                                      'in the same format')
commands.sorttfgenes.outfile.desc = ('The output file. If this is specified '
                                     'with a full path, `--outdir` will be '
                                     'ignored')
commands.sorttfgenes.outfile.callback = (
    lambda opt, ps: opt.set_value(Path(ps.outdir.value).joinpath(opt.value))
    if ps.outdir.value and opt.value and '/' not in opt.value
    else opt.set_value(Path(ps.outdir.value).joinpath(
        Path(ps.origin).with_suffix('.sorted' + Path(ps.origin.value).suffix)
    ))
    if ps.outdir.value and not opt.value
    else None
    if '/' in opt.value
    else 'Either --outdir or --outfile must be specified'
)

### Pair SNP gene using distance restriction
commands.pairsnpgenes = ('Pair SNP-Gene using if SNPs locates in '
                         'flanking regions of genes')
commands.pairsnpgenes.snpbed = GTQC.snpbed
commands.pairsnpgenes.gfile = commands.pwmscan.gfile
commands.pairsnpgenes.inopts = commands.pwmscan.inopts
commands.pairsnpgenes.gcol = commands.pwmscan.gcol
commands.pairsnpgenes.region = commands.pwmscan.region
commands.pairsnpgenes.outfile.desc = ('The output file. If this is specified '
                                      'with a full path, `--outdir` will be '
                                      'ignored')
commands.pairsnpgenes.outfile.callback = (
    lambda opt, ps: opt.set_value(Path(ps.outdir.value).joinpath(opt.value))
    if ps.outdir.value and opt.value and '/' not in opt.value
    else opt.set_value(Path(ps.outdir.value).joinpath(
        Path(ps.snpbed).with_suffix('.snp-gene.txt')
    ))
    if ps.outdir.value and not opt.value
    else None
    if '/' in opt.value # or not ps.outdir.value
    else 'Either --outdir or --outfile must be specified'
)

### Adjust covariates
commands.decov = "Adjust covariates from a set of files with columns variantes and rows samples"
commands.decov.expr.required = True
commands.decov.expr.desc = 'The expression matrix'
commands.decov.covfiles.required = True
commands.decov.covfiles.type = list
commands.decov.covfiles.desc = 'The files with covariates.'
commands.decov.vars.required = True
commands.decov.vars.type = list
commands.decov.vars.desc = 'The variable names to select'
commands.decov.tool = 'linear'
commands.decov.tool.desc = 'The tool used to adjust covariates, linear or peer.'
commands.decov.outfile.required = True
commands.decov.outfile.desc = (
    'The output file. If this is specified with a full path, '
    '`--outdir` will be ignored. Otherwise the outfile will be '
    '`<outdir>/<outfile>`'
)
commands.decov.peerR = params.Rscript.value
commands.decov.peerR.desc = 'The Rscript for peer. Peer is an old package, it only runs on older R'
commands.decov.outfile.required = True
commands.decov.outfile.callback = (
    lambda opt, ps: opt.set_value(Path(ps.outdir.value).joinpath(opt.value))
    if ps.outdir.value and opt.value and '/' not in opt.value
    else None
    if '/' in opt.value or not ps.outdir.value
    else 'Either --outdir or --outfile must be specified'
)

### Call cis-eqtls
commands.eqtl = "Do cis-eQTL calling"
commands.eqtl.gtype.required = True
commands.eqtl.gtype.desc = 'The genotype matrix'
commands.eqtl.expr.required = True
commands.eqtl.expr.desc = 'The expression matrix'
commands.eqtl.genepos = params.refgene
commands.eqtl.genepos.show = True
commands.eqtl.snppos.required = True
commands.eqtl.snppos.desc = 'The SNP coordinates in BED or VCF'
commands.eqtl.dist = 250000
commands.eqtl.dist.desc = 'The distance of a cis-eQTL to its target gene'
commands.eqtl.pval = .05
commands.eqtl.pval.desc = 'The pvalue cutoff'
commands.eqtl.outfile = commands.decov.outfile

### Aggregate results for each SNP
commands.aggr = "Aggregate ceQTL results for each SNP"
commands.aggr.infile.required = True
commands.aggr.infile.desc = 'The ceQTL result file'
commands.aggr.on = "Padj"
commands.aggr.on.desc = "Which column of values to aggregate for the same SNP"
commands.aggr.method = "min"
commands.aggr.method.desc = "The method used to aggregate the pvalues (min or fisher)"
commands.aggr.outfile = commands.decov.outfile

### Q-Q plot for the results
commands.qq = "Q-Q plot for the ceQTL results"
commands.qq.col = "Padj"
commands.qq.col.desc = "Which column of pvalues to use."
commands.qq.outfile = commands.decov.outfile
commands.qq.infile.required = True
commands.qq.infile.desc = 'The ceQTL result file'

### Manhattan plot for the results
commands.man = "Manhattan plot for data"
commands.man.col = "Padj"
commands.man.col.desc = "The column of pvalue or adjusted pvalue to use"
commands.man.infile = commands.qq.infile
commands.man.pos.required = True
commands.man.pos.desc = "The position for each row item in BED"
commands.man.hilight.desc = ["The file used to highlight on the plot",
                             "- It could just one column, then the SNPs will be just "
                             "colored differently",
                             "- Or 2 columns, with SNPs and the labels for each one of them."]
commands.man.outfile = commands.decov.outfile
commands.man.params.type = dict
commands.man.params.desc = "Parameters for `CMplot`. See `?CMplot` in R for details"

### Compare results for different measurements
commands.compare = 'Rank compare of two measurements (i.e. ceQTL, eQTL, median, etc.)'
commands.compare.in1.required = True
commands.compare.in1.desc = "The first input file"
commands.compare.in2.required = True
commands.compare.in2.desc = "The second input file"
commands.compare.col1.required = True
commands.compare.col1.desc = ["The column of 1st file to use for compare (i.e. p/padj values)",
                              "Could be index (0-based) or column name"]
commands.compare.col2.required = True
commands.compare.col2.desc = "The column of 2nd file to use for compare (i.e. p/padj values)"
commands.compare.cut1 = 0.05
commands.compare.cut1.desc = "The cutoff for 1st measurement (will keep records with col1 < cut1)"
commands.compare.cut2 = 0.05
commands.compare.cut2.desc = "The cutoff for 2nd measurement (will keep records with col1 < cut2)"

### Plot ROC using gold standard
commands.roc = 'Plot ROC curves for different measurements'
commands.roc.infiles.type = list
commands.roc.infiles.required = True
commands.roc.infiles.desc = 'The input files'
commands.roc.infiles.callback = (lambda opt, ps: '--outdir is required.'
                                 if not ps.outdir.value else None)
commands.roc.cols.type = list
commands.roc.cols.required = True
commands.roc.cols.desc = 'The column of data for each file used for ranking, usually pvalues'
commands.roc.cols.callback = (lambda opt, ps: "Unequal lengths for infiles and cols"
                              if len(opt.value) != len(ps.infiles.value)
                              else None)
commands.roc.names = []
commands.roc.names.type = list
commands.roc.names.desc = ('Names for each infile to show in the ROC plot. By default, '
                           'it will be stems of input files.')
commands.roc.names.callback = (lambda opt, ps: "Unequal lengths for infiles and names"
                               if opt.value and len(opt.value) != len(ps.infiles.value)
                               else opt.set_value([Path(infile).stem
                                                   for infile in ps.infiles.value])
                               if not opt.value
                               else None)
commands.roc.cuts = []
commands.roc.cuts.type = list
commands.roc.cuts.desc = 'Keep records with `cols` value < `cuts` for each infile'
commands.roc.cuts.callback = (lambda opt, ps: "Unequal lengths for infiles and cuts"
                              if opt.value and len(opt.value) != len(ps.infiles.value)
                              else None)
commands.roc.sep = True
commands.roc.sep.desc = ('Whether draw ROC curves for each infile separately or in one plot. '
                         'If in one plot, it will be save in '
                         '`<outdir>/<names[0]-names[1]-...>.png`'
                         'Otherwise, `<outdir>/<names[0]>.png`, ...')
commands.roc.rev = True
commands.roc.rev.desc = 'Reversely rank the records (after `cuts` applied)'
commands.roc.gold.required = True
commands.roc.gold.desc = 'The gold standard SNPs, one per line.'

### Hypergeometric distrition pvalue for results hitting gold standard
commands.hg = 'Hypergeometric test for measurements hitting gold standard'
commands.hg.infile.required = True
commands.hg.infile.desc = 'The input file'
commands.hg.col = 'Padj'
commands.hg.col.desc = 'The column used to filter records, usually pvalues'
commands.hg.cut = 0.05
commands.hg.cut.desc = 'Only records kept with `col` < `cut`'
commands.hg.bign.required = True
commands.hg.bign.desc = ('The total number of records (background), '
                         'or a file with number of lines equals to that number.')
commands.hg.gold = commands.roc.gold

### Mediation analysis
commands.med = 'Mediation analysis using SNPs as mediator'
commands.med.expr.required = True
commands.med.expr.desc = ('The expression matrix, '
                          'with genes as rows and samples as columns.')
commands.med.gtype.required = True
commands.med.gtype.desc = ('The genotype matrix, '
                     'with SNPs as rows and samples as columns.')
commands.med.tft.required = True
commands.med.tft.desc = ('The TF target matrix, '
                         'with target genes as rows and TFs as columns.')
commands.med.snpgene.required = True
commands.med.snpgene.desc = [
    'The SNP-gene pairs to limit the number of regressions.',
    'This is usually decided by the distance of the SNP to the gene.',
    'You can use `ceqtl-tools snpgene` to generate this file.'
]
commands.med.pcut = 0.05
commands.med.pcut.desc = 'The pvalue cutoff for the trios.'
commands.med.pval = commands.med.pcut
commands.med.njobs = 1
commands.med.njobs.desc = ('Split the cases into different jobs. '
                           'You can distribute the jobs into clusters, for example.')
commands.med.runner = 'local'
commands.med.runner.desc = ('The runner for the pipeline. '
                            'See https://pyppl.readthedocs.io/en/latest/runners/')
commands.med.outfile = commands.decov.outfile

commands._.outdir.desc = 'The output directory'
commands._.nthread = 1
commands._.nthread.desc = 'Number of threads to use'
commands._.ncores = commands._.nthread

def main():
    """Main function"""
    #sys.path.insert(1, str(Path(__file__).parent))
    command, opts, _ = commands._parse(dict_wrapper=Diot)
    modname = command.replace('-', '_')
    module = __import__('ceqtl.tools', fromlist=[modname])
    getattr(module, modname).main(opts)

if __name__ == "__main__":
    main()
