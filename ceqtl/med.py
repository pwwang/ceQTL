import pandas as pd
from pipen import Pipen
from pipen.channel import expand_dir
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.stats import Mediation as Mediation_, MetaPvalue1 as MetaPvalue1_

from shared.procs import MergeChunks


parser.add_argument(
    "--expr",
    help=(
        "The expression matrix file, with genes as columns and samples as rows. "
    ),
    required=True,
    type="path",
)
parser.add_argument(
    "--geno",
    help=(
        "The genotype matrix file, with SNP as columns and samples as rows. "
        "It should have the same number and order of samples as "
        "the expression matrix."
    ),
    required=True,
    type="path",
)
parser.add_argument(
    "--cov",
    help=(
        "The file with covariate variables with covariates as columns and "
        "samples as rows. "
    ),
    default=None,
    type="path",
)
parser.add_argument(
    "--trios",
    help="The ceQTL trio file",
    required=True,
    type="path",
)
parser.add_argument(
    "--nchunks",
    help="Break the trios into chunks to run in parallel",
    default=4,
)


def get_pipeline(args):
    class ComposeInputFiles(Proc):
        """Compose the input files for mediation analysis

        Input:
            exprfile: The expression matrix file
            genofile: The genotype matrix file
            covfile: The file with covariate variables
            triofile: The ceQTL trio file

        Output:
            outdir: The directory to store the input files for mediation analysis

        Envs:
            transpose_expr (flag): Transpose the expression matrix
            transpose_geno (flag): Transpose the genotype matrix
            transpose_cov (flag): Transpose the covariate matrix
        """
        input = "exprfile:file, genofile:file, covfile:file, triofile:file"
        input_data = [(args.expr, args.geno, args.cov, args.trios)]
        output = "outdir:dir:{{in.exprfile | stem0}}.chunks"
        lang = config.lang.rscript
        envs = {
            "transpose_expr": False,
            "transpose_geno": False,
            "transpose_cov": False,
            "nchunks": args.nchunks,
        }
        script = """
            source("{{biopipen_dir}}/utils/misc.R")
            library(dplyr)
            exprfile <- {{in.exprfile | r}}
            genofile <- {{in.genofile | r}}
            covfile <- {{in.covfile | r}}
            triofile <- {{in.triofile | r}}
            outdir <- {{out.outdir | r}}
            transpose_expr <- {{envs.transpose_expr | r}}
            transpose_geno <- {{envs.transpose_geno | r}}
            transpose_cov <- {{envs.transpose_cov | r}}
            nchunks <- {{envs.nchunks | int}}

            expr <- read.table(exprfile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
            geno <- read.table(genofile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
            if (!is.null(covfile)) cov <- read.table(covfile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

            trios <- read.table(triofile, header=TRUE, sep="\t", check.names=FALSE)
            if (transpose_expr) expr <- t(expr)
            if (transpose_geno) geno <- t(geno)
            if (transpose_cov && !is.null(covfile)) cov <- t(cov)

            common_samples = intersect(rownames(expr), rownames(geno))
            if (!is.null(covfile)) common_samples = intersect(common_samples, rownames(cov))

            data <- cbind(expr[common_samples, , drop=FALSE], geno[common_samples, , drop=FALSE])
            if (!is.null(covfile)) data <- cbind(data, cov[common_samples, , drop=FALSE])

            all_trios_idx <- cut(1:nrow(trios), nchunks, labels=FALSE)
            all_trios <- split(trios, all_trios_idx)

            allcovs <- if(is.null(covfile)) NULL else paste(colnames(cov), collapse = ",")
            covs <- if(is.null(covfile)) NULL else colnames(cov)
            for (i in seq_along(all_trios)) {
                log_info("Processing chunk {i} ...")
                datadir <- file.path(outdir, sprintf("chunk-%d", i))
                dir.create(datadir, showWarnings=FALSE, recursive=TRUE)
                chunk_trios <- all_trios[[i]]
                vars <- unique(c(chunk_trios$SNP, chunk_trios$TF, chunk_trios$Target, covs))
                chunk_data <- data[, vars, drop=FALSE]

                datafile <- file.path(datadir, "data.txt")
                write.table(chunk_data, datafile, sep="\t", quote=FALSE, row.names=FALSE)

                casefile <- file.path(datadir, "case.txt")
                casedata <- data.frame(
                    Case = paste0("Case", 1:nrow(chunk_trios)),
                    M = chunk_trios$SNP,
                    X = chunk_trios$TF,
                    Y = chunk_trios$Target,
                    Model_M = "lm",
                    Model_Y = "lm"
                )
                if (!is.null(covfile)) {
                    casedata$Cov <- allcovs
                }
                write.table(casedata, casefile, sep="\t", quote=FALSE, row.names=FALSE)
            }
        """

    class Mediation(Mediation_):
        requires = ComposeInputFiles
        input_data = lambda ch: pd.concat(
            [
                expand_dir(ch, 0, "chunk-*/data.txt"),
                expand_dir(ch, 0, "chunk-*/case.txt")
            ],
            axis=1,
        )

    class MergeMediationChunks(MergeChunks):
        requires = Mediation
        output = "outfile:file:mediation.txt"
        input_data = lambda ch: [ch.outfile.to_list()]
        export = True

    class AdjustPvalue(MetaPvalue1_):
        requires = MergeMediationChunks
        output = "outfile:file:mediation-adjusted.txt"
        envs = {
            "id_cols": ["M", "X", "Y"],
            "pval_col": "Pval",
            "padj": "fdr",
        }

    class MediationAnalysis(Pipen):
        starts = ComposeInputFiles

    return MediationAnalysis()


def main():
    args = parser.parse_extra_args()
    pipen = get_pipeline(args)
    pipen.run()
