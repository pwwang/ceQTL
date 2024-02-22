"""Define processes for the pipeline"""
from pathlib import Path
from datar import f
from datar.base import paste0
from datar.dplyr import mutate, select
from pipen.channel import expand_dir
from pipen.utils import is_loading_pipeline
from pipen_args import parser, config
from pipen_args.parser_ import FallbackNamespace
from biopipen.core.config import config as bp_config
from biopipen.core.proc import Proc
from biopipen.ns.stats import (
    ChowTest as _ChowTest,
    LiquidAssoc as _LiquidAssoc,
    MetaPvalue as _MetaPvalue,
)

from shared.procs import DataPreparation as _DataPreparation
from .args_def import add_args

# Add extra arguments to the parser
parser = add_args(parser)
args = parser.parse_extra_args()
is_loading = is_loading_pipeline() or isinstance(args, FallbackNamespace)


class DataPreparation(_DataPreparation):
    input_data = [(args.geno, args.expr, args.cov, args.snpgene, args.tftarget)]
    envs = {"nchunks": args.nchunks, "ncores": args.ncores}


class MergeChunks(Proc):
    """Merge the chunks of the results.

    Input:
        infiles: The output files of chunks from each model.

    Output:
        outfile: The merged output file.
    """
    input = "infiles:files"
    # outfile = ...  # to be set in subclasses
    lang = bp_config.lang.bash
    script = """
        infile0={{in.infiles | first | quote}}
        outfile={{out.outfile | quote}}
        head -n 1 $infile0 > $outfile
        {% for infile in in.infiles %}
        tail -n +2 {{infile | quote}} >> $outfile
        {% endfor %}
    """


ModelProcs = []
if "ChowTest" in config or is_loading:
    class ChowTest(_ChowTest):
        requires = DataPreparation
        input_data = lambda ch: (
            expand_dir(ch, pattern="chunk-*")
            >> mutate(
                infile=paste0(f.outdir, "/expression.txt"),
                groupfile=paste0(f.outdir, "/genotype.txt"),
                fmlfile=paste0(f.outdir, "/formula-chow.txt"),
            )
            >> select(~f.outdir)
        )
        envs = {"padj": "none"}

    class MergeChowTestChunks(MergeChunks):
        requires = ChowTest
        output = "outfile:file:chowtest.txt"
        input_data = lambda ch: [
            list(
                Path(ch.outfile.iloc[0])
                .parent
                .parent
                .parent
                .glob("*/output/expression.chowtest.txt")
            )
        ]

    ModelProcs.append(ChowTest)

if "InteractionLm" in config or is_loading:
    class InteractionLm(Proc):
        """Run MiRNA-SNP model to find ceQTLs.

        See also https://github.com/gawilk/miRNA-SNP

        It a liner model: target ~ cov + TF * SNP.
        followed by anova to test the significance of the interaction term.

        Input:
            infile: The expression matrix file, including the covariates.
            groupfile: The genotype matrix file.
            fmlfile: The formula file.

        Output:
            outfile: The output file.

        Envs:
            padj (choice): The p-value adjustment method to use.
                - none: No p-value adjustment (no Padj column in outfile).
                - holm: Holm-Bonferroni method.
                - hochberg: Hochberg method.
                - hommel: Hommel method.
                - bonferroni: Bonferroni method.
                - BH: Benjamini-Hochberg method.
                - BY: Benjamini-Yekutieli method.
                - fdr: FDR correction method.
        """
        requires = DataPreparation
        input = "infile:file, groupfile:file, fmlfile:file"
        input_data = lambda ch: (
            expand_dir(ch, pattern="chunk-*")
            >> mutate(
                infile=paste0(f.outdir, "/expression.txt"),
                groupfile=paste0(f.outdir, "/genotype.txt"),
                fmlfile=paste0(f.outdir, "/formula-interactionlm.txt"),
            )
            >> select(~f.outdir)
        )
        output = "outfile:file:{{in.infile | stem}}.interactionlm.txt"
        lang = bp_config.lang.rscript
        envs = {"padj": "none"}
        script = "file://scripts/InteractionLm.R"

    class MergeInteractionLmChunks(MergeChunks):
        requires = InteractionLm
        output = "outfile:file:InteractionLm.txt"
        input_data = lambda ch: [
            list(
                Path(ch.outfile.iloc[0])
                .parent
                .parent
                .parent
                .glob("*/output/expression.interactionlm.txt")
            )
        ]

    ModelProcs.append(InteractionLm)

if "LiquidAssoc" in config or is_loading:
    class LiquidAssoc(_LiquidAssoc):
        requires = DataPreparation
        input_data = lambda ch: (
            expand_dir(ch, pattern="chunk-*")
            >> mutate(
                infile=paste0(f.outdir, "/expression-nocov.txt"),
                covfile=args.cov,
                groupfile=paste0(f.outdir, "/genotype.txt"),
                fmlfile=paste0(f.outdir, "/formula-la.txt"),
            )
            >> select(~f.outdir)
        )
        envs = {"padj": "none", "xyz_names": ["TF", "Target", "SNP"]}

    class MergeLiquidAssocChunks(MergeChunks):
        requires = LiquidAssoc
        output = "outfile:file:liquidassoc.txt"
        input_data = lambda ch: [
            list(
                Path(ch.outfile.iloc[0])
                .parent
                .parent
                .parent
                .glob("*/output/expression-nocov.liquidassoc.txt")
            )
        ]

    ModelProcs.append(LiquidAssoc)

if config and not ModelProcs:
    raise ValueError("No model specified in the config file")


class MetaPvalue(_MetaPvalue):
    requires = ModelProcs
    input_data = lambda *chs: [[ch.iloc[0, 0] for ch in chs]]
    envs = {
        "id_cols": ["SNP", "TF", "Target"],
        "pval_cols": "Pval",
        "padj": "BH",
    }
