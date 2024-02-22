from pipen import Pipen
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.rnaseq import UnitConversion as _UnitConversion


parser.add_argument(
    "--infile",
    help="The input expression file, with genes as rows and samples as columns.",
    required=True,
)
parser.add_argument(
    "--inunit",
    help=(
        "The input unit of the expression values.\n"
        "You can also use an expression to indicate the input unit, e.g., "
        "`log2(counts + 1)`. The expression should be like `A * fn(B*X + C) + D`, "
        "where `A`, `B`, `C` and `D` are constants, `fn` is a function, and X is\n"
        "the input unit.\n"
        "Currently only `expr`, `sqrt`, `log2`, `log10` and `log` are supported as "
        "functions.\n"
        "Supported input units are:\n"
        "* counts/count/rawcounts/rawcount: raw counts.\n"
        "* cpm: counts per million.\n"
        "* fpkm/rpkm: fragments per kilobase of transcript per million.\n"
        "* fpkmuq/rpkmuq: upper quartile normalized FPKM/RPKM.\n"
        "* tpm: transcripts per million.\n"
        "* tmm: trimmed mean of M-values.\n"
        "This is a shotcut for `--UnitConversion.envs.inunit`"
    ),
    required=True,
)
parser.add_argument(
    "--outunit",
    help=(
        "The output unit of the expression values. An expression can also be used for "
        "transformation (e.g. `log2(tpm + 1)`). If `inunit` is `count`, then this "
        "means we are converting raw counts to tpm, and transforming it to "
        "`log2(tpm + 1)` as the output. Any expression supported by `R` can be used. "
        "Same units as `inunit` are supported.\n"
        "This is a shotcut for `--UnitConversion.envs.outunit`"
    ),
    required=True,
)
parser.add_argument(
    "--genefilter",
    help=(
        "An R expression that takes a vector of gene expression values and returns "
        "a logical vector of the same length indicating whether the gene should be "
        "kept or not, after unit conversion. For example, `mean(x) > 1` will keep "
        "genes with mean expression greater than 1.\n"
        "This is a shotcut for `--ExpressionQC.envs.genefilter`"
    ),
    default="TRUE",
)
args = parser.parse_extra_args()


class UnitConversion(_UnitConversion):
    input_data = [args.infile]
    envs = {
        "inunit": args.inunit,
        "outunit": args.outunit,
    }


class ExpressionQC(Proc):
    """Filter genes after unit conversion and transpose the matrix

    Input:
        infile: The input file

    Output:
        outfile: The output file.
            The rows are samples and the columns are genes.

    Envs:
        genefilter: An R expression to filter genes
    """
    requires = UnitConversion
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.rscript
    envs = {"genefilter": args.genefilter}
    script = """
        library(dplyr)
        data <- read.table({{in.infile | r}}, header=TRUE, row.names=1)
        keep <- apply(data, 1, function(x) {{genefilter}})
        data <- data[keep, , drop=FALSE]
        write.table(t(data), {{out.outfile | r}}, sep="\t", quote=FALSE)
    """


class Pipeline(Pipen):
    name = "ceqtl_exprqc"
    desc = __doc__
    starts = UnitConversion


def main():
    Pipeline().run()
