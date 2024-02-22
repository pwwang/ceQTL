from pathlib import Path
from pandas import DataFrame
from pipen import Pipen
from pipen.channel import expand_dir
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.snp import PlinkSimulation
from biopipen.ns.stats import DiffCoexpr as _DiffCoexpr

from shared.procs import DataPreparation


parser.add_argument(
    "-g",
    "--ngenes",
    help="Number of genes to simulate",
    required=True,
    type=int,
)
parser.add_argument(
    "-s",
    "--nsamples",
    help="Number of samples to simulate",
    required=True,
    type=int,
)
parser.add_argument(
    "-p",
    "--nsnps",
    help="Number of SNPs to simulate",
    required=True,
    type=int,
)
parser.add_argument(
    "-n",
    "--ncores",
    help="Number of cores to use for all processes",
    default=config.misc.ncores,
    type=int,
)
parser.add_argument(
    "-c",
    "--nchunks",
    help="Number of chunks to split the data into",
    default=1,
    type=int,
)

args = parser.parse_extra_args(kwargs={"nsamples": 100})


class GenoSimulation(PlinkSimulation):
    input_data = [(args.nsnps, args.nsamples / 2, args.nsamples / 2)]
    envs = {"sample_prefix": "Sample", "seed": 8525}
    export = True


class DiffCoexprSimulation(Proc):
    """Simulate gene expressions with differential co-expression"""
    requires = GenoSimulation
    lang = config.lang.rscript
    input = "genofile:file"
    input_data = lambda ch: ch.gtmat
    output = [
        "exfile:file:simulated_expr.txt",
        "snpgenefile:file:simulated_snpgenes.gmt",
        "tftargetfile:file:simulated_tftargets.gmt"
    ]
    script = "file://scripts/DiffCoexprSimulation.R"
    envs = {
        "seed": 8525,
        "ncores": args.ncores,
        "ngenes": args.ngenes,
        "npair_range": [10, 16],
        "nregulator_range": [.1, .2],
        "ntarget_range": [10, 25],
        "gene_prefix": "Gene",
        "gene_index_start": 1,
        "bionoise": 0.05,
        "expnoise": 0.05,
        "transpose_exprs": True,
    }
    export = True


class ParallelizationPreparation(DataPreparation):
    requires = GenoSimulation, DiffCoexprSimulation
    input_data = lambda ch1, ch2: [
        (
            ch1.gtmat.iloc[0],
            ch2.exfile.iloc[0],
            None,
            ch2.snpgenefile.iloc[0],
            ch2.tftargetfile.iloc[0],
        )
    ]
    lang = config.lang.rscript
    envs = {"ncores": args.ncores, "nchunks": args.nchunks, "transpose_geno": True}


class DiffCoexpr(_DiffCoexpr):
    requires = ParallelizationPreparation
    input_data = lambda ch: DataFrame(
        {
            "infile": expand_dir(ch, pattern="chunk-*").iloc[:, 0].str.cat(
                ["expression.txt"] * args.nchunks,
                sep="/",
            ),
            "groupfile": expand_dir(ch, pattern="chunk-*").iloc[:, 0].str.cat(
                ["genotype.txt"] * args.nchunks,
                sep="/",
            ),
        }
    )
    envs = {
        "ncores": args.ncores,
        "seed": 8525,
        "cor_method": "pearson",
        "beta": 6,
        "padj": "none",
        # "keep": "Padj < 0.05",
        "perm_batch": 20,
    }


class PAdjust(Proc):
    """Concatenate the results and adjust pvalues

    Input:
        infiles: The output files of all chunks.

    Output:
        outfile: The output file with adjusted pvalues.

    Envs:
        padj (choice): The method to adjust pvalues.
            - BH: Benjamini-Hochberg
            - bonferroni: Bonferroni
            - fdr: FDR
            - holm: Holm
            - hochberg: Hochberg
            - hommel: Hommel
            - BY: Benjamini-Yekutieli
    """
    requires = DiffCoexpr
    input = "infiles:files"
    input_data = lambda ch: [
        list(Path(ch.outfile.iloc[0]).parent.parent.glob("*/expression.diffcoexpr.txt"))
    ]
    output = "outfile:file:triplets.txt"
    lang = config.lang.rscript
    envs = {"padj": "BH"}
    script = """
        infiles = {{in.infiles | r}}
        outfile = {{out.outfile | r}}
        method = {{envs.padj | r}}
        out = c()
        for (infile in infiles) {
            indata = read.table(infile, header = TRUE, sep = "\t", row.names = NULL, check.names = FALSE)
            out = rbind(out, indata)
        }
        out$Padj = p.adjust(out$Pval, method)
        out = out[order(out$Padj), ]
        write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    """


class Pipeline(Pipen):
    name = "ceqtl_simulate"
    desc = "Simulate data for ceqtl pipeline"
    starts = GenoSimulation


def main():
    Pipeline().run()
