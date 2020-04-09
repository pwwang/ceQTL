# testing stage
set ceqtl_tools python ceqtl/tool.py
set ceqtl python ceqtl/__init__.py
# for production
# set ceqtl_tools ceqtl-tools
# set ceqtl ceqtl

# -1. download data and preprocess data
#    will have these 4 files:
#    - TCGA-LUAD.survival.tsv.gz
#    - TCGA-LUAD.GDC_phenotype.tsv.gz
#    - TCGA-LUAD.expr.txt # log2(count+1) expression with gene symbols
#    - TCGA-LUAD.gt.txt # genotype matrix with rs ids
$ceqtl_tools rs2bed tcga \
    --cancer LUAD \
    --snpmani TCGA/LUAD/LUAD.gt.manifest.txt \
    --snpmeta TCGA/LUAD/LUAD.gt.meta.json \
    --outdir TCGA/

# 0. Get SNP coordinates
$ceqtl_tools rs2bed \
    --rsfile TCGA/LUAD/TCGA-LUAD.gt.txt \
    --outfile TCGA/LUAD/0.TCGA-LUAD.gt.bed \
    --inopts.skip 1

# 1. Genotype data QC
$ceqtl_tools genotype-qc \
    --gtmat TCGA/LUAD/TCGA-LUAD.gt.txt \
    --outdir TCGA/LUAD/1.genotype-qc \
    --snpbed TCGA/LUAD/0.TCGA-LUAD.gt.bed

# 2. Expression data QC
$ceqtl_tools expr-qc \
    --outdir TCGA/LUAD/2.expr-qc \
    --expr TCGA/LUAD/TCGA-LUAD.expr.txt \
    --intsform 'function(x) 2^x - 1' \
    --inunit count \
    --outunit tmm \
    --premed 1e-9 \
    --postmed 1

# 3. PWM scan the mutated sequence of gene promoters
$ceqtl_tools pwmscan \
    --gfile TCGA/LUAD/2.expr-qc/TCGA-LUAD.expr.txt \
    --outdir TCGA/LUAD/3.pwmscan \
    --inopts.skip 1 \
    --mut TCGA/LUAD/1.genotype-qc/TCGA-LUAD.snp.bed \
    --nthread 60

# 4. Add TF-gene pairs from wildtype sequence motif scan
$ceqtl_tools sorttfgenes \
    --origin TCGA/LUAD/3.pwmscan/TCGA-LUAD.tf-gene.txt \
    --addition ~/reference/motif_databases/wg.250K/tf-genes-q5e-2.txt \
    --outfile TCGA/LUAD/4.TCGA-LUAD.tf-gene.txt

# 5. Pair SNP-gene based on coordinates
$ceqtl_tools pairsnpgenes \
    --snpbed TCGA/LUAD/1.genotype-qc/TCGA-LUAD.snp.bed \
    --gfile TCGA/LUAD/2.expr-qc/TCGA-LUAD.expr.txt \
    --inopts.skip 1 \
    --region 250000 \
    --outfile TCGA/LUAD/5.TCGA-LUAD.snp-gene.txt

# 6. Population structure by PCA the genotype matrix
bioprocs plink.pPlinkPCA -i.indir TCGA/LUAD/1.genotype-qc/TCGA-LUAD.gt.plink/ \
    --args.samid iid --args.select 2 --args.nthread 20 \
    --config.export_dir TCGA/LUAD/6.TCGA-LUAD.popstruct

# 7. Remove covariate effect from expression data
$ceqtl_tools decov \
    --expr TCGA/LUAD/2.expr-qc/TCGA-LUAD.expr.txt \
    --covfiles TCGA/LUAD/TCGA-LUAD.GDC_phenotype.tsv.gz \
               TCGA/LUAD/6.TCGA-LUAD.popstruct/TCGA-LUAD.gt.pcs.txt \
    --outfile TCGA/LUAD/7.expr.decoved.txt \
    --vars ROWNAME PC1 PC2 submitter_id.samples days_to_birth.demographic \
           ethnicity.demographic gender.demographic tumor_stage.diagnoses

$ceqtl_tools decov \
    --expr TCGA/LUAD/2.expr-qc/TCGA-LUAD.expr.txt \
    --covfiles TCGA/LUAD/TCGA-LUAD.GDC_phenotype.tsv.gz \
               TCGA/LUAD/6.TCGA-LUAD.popstruct/TCGA-LUAD.gt.pcs.txt \
    --outfile TCGA/LUAD/7.expr.peer.txt \
    --vars ROWNAME PC1 PC2 submitter_id.samples days_to_birth.demographic \
           ethnicity.demographic gender.demographic tumor_stage.diagnoses \
    --tool peer
    --peerR ~/miniconda3/envs/r-peer/bin/Rscript

# 8. Run ceqtl
$ceqtl --expr TCGA/LUAD/7.expr.decoved.txt \
    --gtype TCGA/LUAD/1.genotype-qc/TCGA-LUAD.gt.txt \
    --tft TCGA/LUAD/4.TCGA-LUAD.tf-gene.txt \
    --njobs 100 \
    --pval 1.1 \
    --snpgene TCGA/LUAD/5.TCGA-LUAD.snp-gene.txt \
    --runner sge1d-thr1 \
    --outfile TCGA/LUAD/8.decov.chow.txt
$ceqtl --expr TCGA/LUAD/7.expr.peer.txt  \
    --gtype TCGA/LUAD/1.genotype-qc/TCGA-LUAD.gt.txt  \
    --tft TCGA/LUAD/4.TCGA-LUAD.tf-gene.txt  \
    --njobs 100 \
    --pval 1.1 \
    --snpgene TCGA/LUAD/5.TCGA-LUAD.snp-gene.txt \
    --runner sge1d-thr1 \
    --outfile TCGA/LUAD/8.peer.chow.txt

# 9. Aggregate pvalues for the same SNPs
$ceqtl_tools aggr --infile TCGA/LUAD/8.decov.chow.txt --outfile 9.decov.aggrchow.txt

# 10. QQ plot
$ceqtl_tools qq --infile TCGA/LUAD/9.decov.aggrchow.txt --outfile TCGA/LUAD/10.qq.png

# 11. Manhattan plot
$ceqtl_tools man --infile TCGA/LUAD/9.decov.aggrchow.txt \
    --pos TCGA/LUAD/1.genotype-qc/TCGA-LUAD.snp.bed \
    --outfile TCGA/LUAD/11.manhattan.jpg \
    --params.ylim:list 1 10

# 12. call eqtl
$ceqtl_tools eqtl \
    --gtype TCGA/LUAD/1.genotype-qc/TCGA-LUAD.gt.txt \
    --expr TCGA/LUAD/7.expr.decoved.txt \
    --snppos TCGA/LUAD/1.genotype-qc/TCGA-LUAD.snp.bed \
    --pval 1 \
    --outfile TCGA/LUAD/12.ciseqtl.txt

# 13. aggreate eqtl calls on each SNP
$ceqtl_tools aggr \
    --infile TCGA/LUAD/12.ciseqtl.txt \
    --outfile TCGA/LUAD/13.ciseqtl.aggr.txt \
    --on p-value

# 14. use eqtl as gold standard to evaluate ceQTL
for pcut in 0.05 0.01 0.005 0.001
    $ceqtl_tools roc \
        --infiles TCGA/LUAD/9.decov.aggrchow.txt \
        --cols Padj \
        --gold (awk "\$4<$pcut" TCGA/LUAD/13.ciseqtl.aggr.txt | cut -f1 | psub) \
        --outdir TCGA/LUAD/14.eqtl.roc-$pcut
end

# 15. Get gwas snps from gwasdb/gwascatlog
cd /path/to/gwasdb
python query.py query "lung cancer" > /path/to/TCGA/LUAD/15.LUAD.gwas.txt
python query.py query "lung adenocarcinoma" | tail -n +2 >> /path/to/TCGA/LUAD/15.LUAD.gwas.txt
grep -i 'lung cancer\|lung carcinoma' gwas_catalog_v1.0.2-associations_e98_r2020-03-08.tsv \
    > /path/to/TCGA/LUAD/15.LUAD.gwascatlog.txt
# get available snps (overlapping with TCGA probes)
cd /path/to/ceQTL
begin
    cut -f3 TCGA/LUAD/15.LUAD.gwas.txt
    cut -f22 TCGA/LUAD/15.LUAD.gwascatlog.txt
end | sort -u | grep '^rs' > TCGA/LUAD/15.LUAD.gwas_snps.txt

# 16. compare recovery of gwas snps
for pcut in 0.05 0.01 0.005 0.001
    $ceqtl_tools roc \
        --infiles TCGA/LUAD/9.decov.aggrchow.txt TCGA/LUAD/13.ciseqtl.aggr.txt \
        --cols Pval p-value \
        --gold TCGA/LUAD/15.LUAD.gwas_snps.txt  \
        --nthread 2 \
        --cuts $pcut $pcut \
        --outdir TCGA/LUAD/16.c_eqtl.roc-$pcut
end

# 17. mediation/moderation analysis
#   mediation
$ceqtl_tools med \
    --expr TCGA/LUAD/2.expr-qc/TCGA-LUAD.expr.txt \
    --gtype TCGA/LUAD/1.genotype-qc/TCGA-LUAD.gt.txt \
    --tft TCGA/LUAD/4.TCGA-LUAD.tf-gene.txt \
    --snpgene TCGA/LUAD/5.TCGA-LUAD.snp-gene.txt \
    --outfile TCGA/LUAD/17.medation.txt \
    --njobs 100 \
    --runner sge1d-thr1
#   moderation
$ceqtl_tools med \
    --expr TCGA/LUAD/2.expr-qc/TCGA-LUAD.expr.txt \
    --gtype TCGA/LUAD/1.genotype-qc/TCGA-LUAD.gt.txt \
    --tft TCGA/LUAD/4.TCGA-LUAD.tf-gene.txt \
    --snpgene TCGA/LUAD/5.TCGA-LUAD.snp-gene.txt \
    --outfile TCGA/LUAD/17.moderation.txt \
    --njobs 100 \
    --runner sge1d-thr1 \
    --type mod

# 18. roc using mediation/moderation as gold standards
$ceqtl_tools roc \
    --infiles TCGA/LUAD/9.decov.aggrchow.txt \
    --cols Padj \
    --gold (cat TCGA/LUAD/17.medation.txt   | cut -f2 | sort -u | psub) \
    --outdir TCGA/LUAD/18.mediation.roc
$ceqtl_tools roc \
    --infiles TCGA/LUAD/9.decov.aggrchow.txt \
    --cols Padj \
    --gold (cat TCGA/LUAD/17.moderation.txt  | cut -f2 | sort -u | psub) \
    --outdir TCGA/LUAD/18.moderation.roc

# 19. atsnp analysis
$ceqtl_tools atsnp \
    --infile TCGA/LUAD/8.decov.chow.txt \
    --snpbed TCGA/LUAD/1.genotype-qc/TCGA-LUAD.snp.bed \
    --outfile TCGA/LUAD/19.atsnp.txt \
    --nthread 20

# 20. use atsnp as gold standard
for pcut in 0.05 0.01 0.005 0.001
    $ceqtl_tools roc \
        --infiles (awk '{print $1"_"$2"\t"$0}' TCGA/LUAD/8.decov.chow.txt | psub) \
        --cols Padj \
        --gold (awk "\$6 < $pcut {print \$1\"_\"\$2\"\t\"\$0}" TCGA/LUAD/19.atsnp.txt | psub) \
        --outdir TCGA/LUAD/20.atsnp.roc-$pcut
end

# 21. pathway enrichment analysis
