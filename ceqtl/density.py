from pipen import Pipen
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import Manhattan as Manhattan_, DensityPlot as DensityPlot_


parser.add_extra_argument(
    "--snp",
    help="The SNP to plot",
)
parser.add_extra_argument(
    "--tf",
    help="The TF to plot",
)
parser.add_extra_argument(
    "--target",
    help="The target gene to plot",
)
parser.add_extra_argument(
    "--expr",
    help=("The expression matrix file, with genes as columns and samples as rows. "),
    type="path",
)
parser.add_extra_argument(
    "--geno",
    help=(
        "The genotype matrix file, with SNP as columns and samples as rows. "
        "It should have the same number and order of samples as "
        "the expression matrix."
    ),
    type="path",
)
parser.add_extra_argument(
    "--cov",
    help=(
        "The file with covariate variables with covariates as columns and "
        "samples as rows. When provided, the target will be residualized "
        "against the covariates first."
    ),
    default=None,
    type="path",
)
parser.add_extra_argument(
    "--transpose-expr",
    dest="transpose_expr",
    help="Transpose the expression matrix",
    action="store_true",
)
parser.add_extra_argument(
    "--transpose-geno",
    dest="transpose_geno",
    help="Transpose the genotype matrix",
    action="store_true",
)
parser.add_extra_argument(
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
        output = "outfile:file:dist-data.txt"
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

            expr <- read.table(exprfile, header=TRUE, row.names=NULL, sep="\\t", check.names=FALSE)
            gene_col <- colnames(expr)[1]
            genes <- make.unique(expr[[gene_col]])
            expr <- expr[, -1, drop=FALSE]
            rownames(expr) <- genes

            if (transpose_expr) {
                expr <- t(expr)
            }
            expr <- expr[, c(tf, target), drop=FALSE]

            geno <- read.table(genofile, header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
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

        class PrepareDataResidualized(PrepareData): ...

    class CalculateCorrs(Proc):
        """Calculate correlations between TF and target by different genotypes

        Input:
            infile: The prepared data file

        Output:
            outfile: The correlations file
        """

        input = "infile:file"
        output = f"outfile:file:{args.snp}-{args.tf}-{args.target}.correlations.txt"
        envs = {
            "snp": args.snp,
            "tf": args.tf,
            "target": args.target,
            "seed": 8525,
            "n_repeats": 100,
            "n_samples_frac": 0.5,
        }
        requires = PrepareData
        lang = config.lang.rscript
        script = """
            infile <- {{in.infile | r}}
            outfile <- {{out.outfile | r}}
            snp <- {{envs.snp | r}}
            tf <- {{envs.tf | r}}
            target <- {{envs.target | r}}
            n_repeats <- {{envs.n_repeats | r}}
            n_samples_frac <- {{envs.n_samples_frac | r}}
            set.seed({{envs.seed | r}})

            data <- read.table(infile, header=TRUE, sep="\\t", check.names=FALSE)
            data <- data[, c(tf, target, snp)]
            # we want to check the distribution of correlations of TF-target in each
            # genotype group. However, there should be only one correlation value
            # per genotype group. So we will sample a subset (half?) of data points in each
            # genotype group to calculate correlations multiple times.
            gts <- unique(data[[snp]])
            out_df <- NULL
            for (gt in gts) {
                n_samples <- ceiling(sum(data[[snp]] == gt) * n_samples_frac)
                corrs <- c()
                for (i in 1:n_repeats) {
                    subset_data <- data[data[[snp]] == gt, , drop=FALSE]
                    sampled_data <- subset_data[sample(nrow(subset_data), n_samples), , drop=FALSE]
                    cor_value <- cor(sampled_data[[tf]], sampled_data[[target]], method="pearson")
                    corrs <- c(corrs, cor_value)
                }
                temp_df <- data.frame(
                    genotype = rep(gt, n_repeats),
                    correlation = corrs
                )
                out_df <- rbind(out_df, temp_df)
            }
            write.table(out_df, outfile, sep="\\t", quote=FALSE, row.names=FALSE)
        """

    if args.cov:

        class CalculateCorrsResidualized(CalculateCorrs):
            requires = PrepareDataResidualized

    class DensityPlot(DensityPlot_):
        envs = {
            "val_col": "correlation",
            "group_by": "genotype",
            "title": f"{args.snp} - {args.tf} - {args.target}",
            "devpars": {"res": 100, "width": 700, "height": 600},
        }
        requires = CalculateCorrs

    if args.cov:

        class DensityPlotResidualized(DensityPlot):
            requires = CalculateCorrsResidualized

    class Pipeline(Pipen):
        starts = PrepareDataResidualized if args.cov else PrepareData
        name = f"Density-{args.snp}-{args.tf}-{args.target}"

    return Pipeline(plugins=["-report"])


def main():
    args = parser.parse_extra_args()
    try:
        pipen = get_pipeline(args)
    except AttributeError as e:
        print("Error in arguments: {}".format(e))
        parser.print_help()
    else:
        pipen.run()
