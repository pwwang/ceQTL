import random
from simpleconf import Config
# from datar import f
from datar.tibble import tibble
# from datar.dplyr import group_by, summarise
from pipen import Pipen
from pipen.channel import expand_dir
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.snp import PlinkSimulation
# from biopipen.ns.stats import DiffCoexpr as _DiffCoexpr

# from shared.procs import DataPreparation


parser.add_argument(
    "--config",
    help=(
        "The configuration file for simulations.\n"
        "The configuration file should be a toml file with the following keys:\n"
        "- seed: The seed for random number generation\n"
        "- nsims: The number of simulations to run\n"
        "- nsnps: The number of SNPs to simulate. 2 numbers to indicate the range.\n"
        "- nsamples: The number of samples to simulate. 2 numbers to indicate the range.\n"
        "- label: The prefix label for SNPs\n"
        "- prevalence: The prevalence of the disease in SNP simulation, 2 numbers to indicate the range.\n"
        "- minfreq: The minimum frequency of the SNP, 2 numbers for range\n"
        "- maxfreq: The maximum frequency of the SNP, 2 numbers for range\n"
        "- hetodds: The odds ratio for heterozygous, 2 numbers for range\n"
        "- homodds: The odds ratio for homozygous, 2 numbers for range\n"
        "- missing: The missing rate of the SNP, 2 numbers for range\n"
        "- sample_prefix: The prefix for the sample names\n"
        "- ngenes: The number of genes to simulate\n"
        "- npair_range: The range of number of gene pairs to simulate\n"
        "- nregulator_range: The range of number of regulators for each gene\n"
        "- ntarget_range: The range of number of targets for each regulator\n"
        "- gene_prefix: The prefix for the gene names\n"
        "- noise: The noise level for gene expression\n"
        "Required keys: nsims, nsnps, nsamples, and ngenes"
    ),
    required=True,
)
parser.add_argument(
    "-c",
    "--nchunks",
    help="Number of chunks to split the data into",
    default=1,
    type=int,
)

args = parser.parse_extra_args(
    kwargs={"config": {"nsims": 10, "nsnps": 100, "nsamples": 100, "ngenes": 100}}
)
conf = Config.load(args.config)
if (
    "nsims" not in conf
    or "nsnps" not in conf
    or "nsamples" not in conf
    or "ngenes" not in conf
):
    raise ValueError(
        "nsims, nsnps, nsamples, and ngenes are required in the configuration file."
    )

seed = conf.get("seed", 8525)
random.seed(seed)
# generate seeds for simulations
seeds = ",".join(str(s) for s in random.sample(range(1, 10000), conf.nsims))


class GenoSimConfigs(Proc):
    """Generate the configurations for the simulations"""
    input = "infile:file, seeds"
    input_data = [(args.config, seeds)]
    output = "outdir:dir:configs"
    lang = config.lang.rscript
    script = """
        library(RcppTOML)
        library(rlang)
        infile = {{in.infile | r}}
        seeds = {{in.seeds | split: "," | r}}
        outdir = {{out.outdir | r}}

        ensure_range = function(x) if(length(x) == 1) c(x, x) else x

        conf = parseTOML(infile)
        nsamples = ensure_range(conf$nsamples)
        nsnps = ensure_range(conf$nsnps)
        nsamples = sample(nsamples[1]:nsamples[2], 1)
        nsnps = sample(nsnps[1]:nsnps[2], 1)

        if (conf$nsamples %% 2 != 0) {
            ncases = floor(conf$nsamples / 2)
            nctrls = ceiling(conf$nsamples / 2)
        } else {
            ncases = conf$nsamples / 2
            nctrls = conf$nsamples / 2
        }
        for (i in 1:conf$nsims) {
            lines = c()
            key = paste0("sim", i)
            lines = c(lines, paste0("nsnps = ", nsnps))
            lines = c(lines, paste0("ncases = ", ncases))
            lines = c(lines, paste0("nctrls = ", nctrls))
            lines = c(lines, paste0("seed = ", seeds[i]))

            if (!is.null(conf$label))
                lines = c(lines, paste0("label = ", conf$label))
            if (!is.null(conf$prevalence)) {
                prevalence = ensure_range(conf$prevalence)
                prevalence = sample(prevalence[1]:prevalence[2], 1)
                lines = c(lines, paste0("prevalence = ", prevalence))
            }
            if (!is.null(conf$minfreq)) {
                minfreq = ensure_range(conf$minfreq)
                minfreq = sample(minfreq[1]:minfreq[2], 1)
                lines = c(lines, paste0("minfreq = ", minfreq))
            }
            if (!is.null(conf$maxfreq)) {
                maxfreq = ensure_range(conf$maxfreq)
                maxfreq = sample(maxfreq[1]:maxfreq[2], 1)
                lines = c(lines, paste0("maxfreq = ", maxfreq))
            }
            if (!is.null(conf$hetodds)) {
                hetodds = ensure_range(conf$hetodds)
                hetodds = sample(hetodds[1]:hetodds[2], 1)
                lines = c(lines, paste0("hetodds = ", hetodds))
            }
            if (!is.null(conf$homodds)) {
                homodds = ensure_range(conf$homodds)
                homodds = sample(homodds[1]:homodds[2], 1)
                lines = c(lines, paste0("homodds = ", homodds))
            }
            if (!is.null(conf$missing)) {
                missing = ensure_range(conf$missing)
                missing = sample(missing[1]:missing[2], 1)
                lines = c(lines, paste0("missing = ", missing))
            }
            if (!is.null(conf$sample_prefix))
                lines = c(lines, paste0("sample_prefix = ", conf$sample_prefix))
            simfile = file.path(outdir, paste0(key, ".toml"))
            writeLines(lines, simfile)
        }
    """


class GenoSimulation(PlinkSimulation):
    requires = GenoSimConfigs
    input_data = lambda ch: expand_dir(ch, pattern="*.toml").iloc[:, 0]
    envs = {
        "sample_prefix": conf.get("seed", "Sample"),
        "seed": conf.get("seed", 8525),
    }
    export = True


class DiffCoexprSimulation(Proc):
    """Simulate gene expressions with differential co-expression

    Input:
        genofile: The genotype file

    Output:
        exfile: The simulated gene expression file
        snpgenefile: The SNP-gene pairs file
        tftargetfile: The TF-target pairs file
        truthallfile: The truth file for all the simulated data
            All trios that are intended to be simulated
        truthfile: The truth file for the simulated data.
            Actual trios with `AnySig` True
            (any of the observed p-values < 0.05 for genotype groups) and
            `MaxObCorrDiff` >= 0.5 (maximum observed correlation difference
            between genotype groups)

    Envs:
        ngenes (type=int): The number of genes to simulate
        npair_range (list): The range of number of gene pairs to simulate
            for each snp.
        nregulator_range (list): The range of number of regulators as proportion
            in the pool.
        ntarget_range (list): The range of number of targets for each regulator.
        gene_prefix: The prefix for the gene names
        gene_index_start (type=int): The starting index for the gene names
        noise (type=float): The noise level for gene expression
        transpose_exprs (flag): Whether to transpose the expression matrix.
    """
    requires = GenoSimulation
    lang = config.lang.rscript
    input = "genofile:file,seed"
    input_data = lambda ch: tibble(genofile=ch.gtmat, seed=seeds.split(","))
    output = [
        "exfile:file:simulated_expr.txt",
        "snpgenefile:file:simulated_snpgenes.gmt",
        "tftargetfile:file:simulated_tftargets.gmt",
        "truthallfile:file:simulated_truth_all.txt",
        "truthfile:file:simulated_truth.txt",
    ]
    script = "file://scripts/DiffCoexprSimulation.R"
    envs = {
        "ncores": 1,
        "ngenes": conf.ngenes,
        "npair_range": conf.get("npair_range", [10, 16]),
        "nregulator_range": conf.get("nregulator_range", [.1, .2]),
        "ntarget_range": conf.get("ntarget_range", [10, 25]),
        "gene_prefix": conf.get("gene_prefix", "Gene"),
        "gene_index_start": 1,
        "noise": conf.get("noise", 0.05),
        "transpose_exprs": True,
    }
    # export = True


# class ParallelizationPreparation(DataPreparation):
#     requires = GenoSimulation, DiffCoexprSimulation
#     input_data = lambda ch1, ch2: tibble(
#         geno=ch1.gtmat,
#         expr=ch2.exfile,
#         cov=None,
#         snpgene=ch2.snpgenefile,
#         tftarget=ch2.tftargetfile,
#     )
#     lang = config.lang.rscript
#     envs = {
#         "ncores": 1,
#         "nchunks": args.nchunks,
#         "transpose_geno": True,
#         "transpose_expr": True,
#     }


# class DiffCoexpr(_DiffCoexpr):
#     requires = ParallelizationPreparation
#     # output = "outfile:file:{{in.infile | stem}}.diffcoexpr.txt"
#     # bring in the dirname of the input file
#     output = [
#         "outfile:file:{{in.infile | dirname | basename}}.diffcoexpr.txt",
#         "group:var:{{in.infile | dirname | stem | replace: '-gtmat', ''}}",
#     ]
#     input_data = lambda ch: tibble(
#         infile=sum(
#             (
#                 expand_dir(tibble(outdir), pattern="*.chunk-*").iloc[:, 0].str.cat(
#                     ["expression.txt"] * args.nchunks,
#                     sep="/",
#                 ).tolist() for outdir in ch.outdir
#             ),
#             [],
#         ),
#         groupfile=sum(
#             (
#                 expand_dir(tibble(outdir), pattern="*.chunk-*").iloc[:, 0].str.cat(
#                     ["genotype.txt"] * args.nchunks,
#                     sep="/",
#                 ).tolist() for outdir in ch.outdir
#             ),
#             [],
#         ),
#     )
#     envs = {
#         "ncores": 1,
#         "seed": 8525,
#         "cor_method": "pearson",
#         "beta": 6,
#         "padj": "none",
#         # "keep": "Padj < 0.05",
#         "perm_batch": 20,
#     }


# class FilterTrios(Proc):
#     """Filter the trios with certain pairs provided by GMTs"""
#     requires = DiffCoexpr, DiffCoexprSimulation
#     input = "infiles:files, snpgene:file, tftarget:file"
#     input_data = lambda ch1, ch2: tibble(
#         ch1
#         >> group_by(f.group)
#         >> summarise(infiles=f.outfile.agg(list)),
#         snpgene=ch2.snpgenefile,
#         tftarget=ch2.tftargetfile,
#     )
#     output = "outfile:file:{{in.snpgene | stem}}.trios.txt"
#     lang = config.lang.rscript
#     script = """
#         infiles = {{in.infiles | r}}
#         sgfile = {{in.snpgene | r}}
#         ttfile = {{in.tftarget | r}}
#         outfile = {{out.outfile | r}}

#         readGMT <- function(gmtfile) {
#             gmt <- list()
#             con <- file(gmtfile, "r")
#             while (TRUE) {
#                 line <- readLines(con, n = 1)
#                 if (length(line) == 0) {
#                     break
#                 }
#                 items <- strsplit(trimws(line), "\t")[[1]]
#                 gmt[[items[1]]] <- items[3:length(items)]
#             }
#             close(con)
#             return(gmt)
#         }

#         df = c()
#         for (infile in infiles) {
#             indata = read.table(
#                 infile,
#                 header = TRUE,
#                 sep = "\t",
#                 row.names = NULL,
#                 check.names = FALSE)
#             df = rbind(df, indata)
#         }
#         df2 = df
#         df2$Feature1 = df$Feature2
#         df2$Feature2 = df$Feature1
#         df = rbind(df, df2)

#         sg = readGMT(sgfile)
#         tt = readGMT(ttfile)

#         # keep only when Feature1 is a TF
#         df = df[df$Feature1 %in% names(tt), ]
#         # keep rows only when Feature2 is a gene of a SNP, where SNP is the name of sg
#         # and gene is one of sg[[SNP]]
#         df$.keep = apply(df, 1, function(x) {
#             x[3] %in% sg[[x[1]]] & x[3] %in% tt[[x[2]]]
#         })
#         df = df[df$.keep, ]
#         df$.keep = NULL

#         write.table(df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
#     """


# class PAdjust(Proc):
#     """Concatenate the results and adjust pvalues

#     Input:
#         infile: The filtered trio file

#     Output:
#         outfile: The output file with adjusted pvalues.

#     Envs:
#         padj (choice): The method to adjust pvalues.
#             - BH: Benjamini-Hochberg
#             - bonferroni: Bonferroni
#             - fdr: FDR
#             - holm: Holm
#             - hochberg: Hochberg
#             - hommel: Hommel
#             - BY: Benjamini-Yekutieli
#     """
#     requires = FilterTrios
#     input = "infile:file"
#     output = "outfile:file:{{in.infile | stem0}}.trios.txt"
#     lang = config.lang.rscript
#     envs = {"padj": "BH"}
#     script = """
#         infile = {{in.infile | r}}
#         outfile = {{out.outfile | r}}
#         method = {{envs.padj | r}}

#         out = read.table(infile, header = TRUE, sep = "\t", row.names = NULL)
#         out$Padj = p.adjust(out$Pval, method)
#         out = out[order(out$Padj), ]
#         write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
#     """


class Pipeline(Pipen):
    name = "ceqtl_simulate"
    desc = "Simulate data for ceqtl pipeline"
    starts = GenoSimConfigs


def main():
    Pipeline().run()
