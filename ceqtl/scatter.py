from pipen import Pipen
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import Manhattan as Manhattan_, Scatter as Scatter_


parser.add_argument(
    "--snp",
    help="The SNP to plot",
    required=True,
)
parser.add_argument(
    "--tf",
    help="The TF to plot",
    required=True,
)
parser.add_argument(
    "--target",
    help="The target gene to plot",
    required=True,
)
parser.add_argument(
    "--split",
    help="Split the plot by the SNP genotypes",
    action="store_true",
)
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
        "samples as rows. When provided, the target will be residualized "
        "against the covariates first."
    ),
    default=None,
    type="path",
)
parser.add_argument(
    "--transpose-expr",
    dest="transpose_expr",
    help="Transpose the expression matrix",
    action="store_true",
)
parser.add_argument(
    "--transpose-geno",
    dest="transpose_geno",
    help="Transpose the genotype matrix",
    action="store_true",
)
parser.add_argument(
    "--transpose-cov",
    dest="transpose_cov",
    help="Transpose the covariate matrix",
    action="store_true",
)


def get_pipeline(args):

    class PrepareData(Proc):
        """Prepare the data for the plot

        Input:
            exprfile: The expression matrix file
            genofile: The genotype matrix file
            covfile: The covariate matrix file
            snp: The SNP to plot
            tf: The TF to plot
            target: The target gene to plot

        Output:
            outfile: The prepared data file
        """
        input = [
            "exprfile:file",
            "genofile:file",
            "covfile:file",
            "snp:var",
            "tf:var",
            "target:var",
        ]
        output = "outfile:file:scatter-data.txt"
        input_data = [(args.expr, args.geno, args.cov, args.snp, args.tf, args.target)]
        lang = config.lang.rscript
        envs = {
            "transpose_expr": args.transpose_expr,
            "transpose_geno": args.transpose_geno,
            "transpose_cov": args.transpose_cov,
        }
        script = """
            library(dplyr)

            exprfile <- {{in.exprfile | r}}
            genofile <- {{in.genofile | r}}
            covfile <- {{in.covfile | r}}
            snp <- {{in.snp | r}}
            tf <- {{in.tf | r}}
            target <- {{in.target | r}}
            outfile <- {{out.outfile | r}}
            transpose_expr <- {{envs.transpose_expr | r}}
            transpose_geno <- {{envs.transpose_geno | r}}
            transpose_cov <- {{envs.transpose_cov | r}}

            expr <- read.table(exprfile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
            if (transpose_expr) {
                expr <- t(expr)
            }
            expr <- expr[, c(tf, target), drop=FALSE]

            geno <- read.table(genofile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
            if (transpose_geno) {
                geno <- t(geno)
            }
            geno <- as.data.frame(geno[, snp, drop=FALSE])
            geno[[snp]] <- case_when(
                geno[[snp]] == 0 ~ "AA",
                geno[[snp]] == 1 ~ "AB",
                geno[[snp]] == 2 ~ "BB",
                TRUE ~ "NA"
            )

            if (!is.null(covfile)) {
                cov <- read.table(covfile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
                if (transpose_cov) {
                    cov <- t(cov)
                }
                common_samples <- intersect(rownames(expr), rownames(geno))
                common_samples <- intersect(common_samples, rownames(cov))
                df <- cbind(expr[common_samples, target, drop=FALSE], cov[common_samples, , drop=FALSE])
                df <- as.data.frame(df)
                fit <- lm(as.formula(sprintf("`%s` ~ .", target)), data=df)
                df[[target]] <- residuals(fit)
                df[[tf]] <- expr[common_samples, tf]
                df <- cbind(df, geno[common_samples, , drop=FALSE])
            } else {
                common_samples <- intersect(rownames(expr), rownames(geno))
                df <- cbind(expr[common_samples, , drop=FALSE], geno[common_samples, , drop=FALSE])
            }

            write.table(df, outfile, sep="\t", quote=FALSE, row.names=FALSE)
        """

    if args.cov:
        class PrepareDataResidualized(PrepareData):
            ...

    class Scatter(Scatter_):
        envs = {"x_col": args.tf, "y_col": args.target}
        envs["stats"] = {
            "poly_line": {},
            "poly_eq": {
                "mapping": """
                    aes(
                        label = paste(
                            after_stat(eq.label),
                            after_stat(rr.label),
                            sep = '*", "*'
                        )
                    )
                """
            },
        }

        if args.cov:
            envs["ggs"] = [f"ylab('{args.target} (residualized)')"]
        else:
            envs["ggs"] = []

        if args.split:
            envs["ggs"].append(f"facet_wrap(~`{args.snp}`)")
        else:
            envs["mapping"] = f"color = `{args.snp}`"

    if args.split and args.cov:
        class ScatterSplitResidualized(Scatter):
            requires = PrepareDataResidualized
    elif args.split:
        class ScatterSplit(Scatter):
            requires = PrepareData
    elif args.cov:
        class ScatterResidualized(Scatter):
            requires = PrepareDataResidualized
    else:
        class Scatter(Scatter):
            requires = PrepareData

    class Pipeline(Pipen):
        starts = PrepareDataResidualized if args.cov else PrepareData
        name = f"Scatter-{args.snp}-{args.tf}-{args.target}"
        outdir = f"Scatter/{args.snp}-{args.tf}-{args.target}"

    return Pipeline()


def main():
    args = parser.parse_extra_args()
    pipen = get_pipeline(args)
    pipen.run()
