from pipen import Pipen
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import Manhattan as Manhattan_


parser.add_extra_argument(
    "--ceqtl-trios",
    dest="ceqtl_trios",
    help="The ceQTL trio file",
    required=True,
)
parser.add_extra_argument(
    "--ceqtl-vars",
    dest="ceqtl_vars",
    help="The ceQTL variant file",
)
parser.add_extra_argument(
    "--bedfile",
    help=(
        "The BED file giving the coordinates of the SNPs. "
        "The name must match the SNP names in the ceQTL file."
    ),
    required=True,
)
parser.add_extra_argument(
    "--chroms",
    help="Order of chromosomes in the plot",
    default="chr1-22",
)
parser.add_extra_argument(
    "--signif",
    help="Significance levels",
    default="5e-8,1e-5",
)
parser.add_extra_argument(
    "--trio-pval-col",
    dest="trio_pval_col",
    help="The column name of the p-values in the ceQTL trio file",
    default="Pval",
)
parser.add_extra_argument(
    "--var-pval-col",
    dest="var_pval_col",
    help="The column name of the p-values in the ceQTL variant file",
    default="MetaPval",
)


def get_pipeline(args):
    class AddCoordsToTrios(Proc):
        """Add coordinates to the ceQTL file"""
        input = "cefile:file, bedfile:file"
        output = "outfile:file:{{in.cefile | basename}}"
        input_data = [(args.ceqtl_trios, args.bedfile)]
        lang = config.lang.rscript
        script = """
            library(dplyr)

            cefile <- {{in.cefile | r}}
            bedfile <- {{in.bedfile | r}}
            outfile <- {{out.outfile | r}}

            cedata <- read.table(cefile, header=TRUE, sep="\\t", check.names=FALSE)
            beddata <- read.table(bedfile, header=FALSE, sep="\\t", check.names=FALSE)
            beddata <- beddata[, 1:4, drop=FALSE]
            colnames(beddata) <- c("chr", "start", "pos", "SNP")
            cedata <- cedata %>% left_join(beddata, by="SNP")

            write.table(cedata, outfile, sep="\\t", quote=FALSE, row.names=FALSE)
        """

    class AddCoordsToVariants(AddCoordsToTrios):
        input_data = [(args.ceqtl_vars, args.bedfile)]


    class ManhattanTrios(Manhattan_):
        requires = AddCoordsToTrios
        envs = {
            "chrom_col": "chr",
            "pos_col": "pos",
            "pval_col": args.trio_pval_col,
            # "label_col": "SNP",
            "chroms": args.chroms,
            "zoom": args.chroms,
            "signif": args.signif,
            # "zoom": ["chr5", "chr7", "chr9", "chr12", "chr15"],
            # "title": f"Manhattan Plot ({args.pval_col})",
            "ylabel": f"-log10({args.trio_pval_col})",
        }

    class ManhattanVariants(Manhattan_):
        requires = AddCoordsToVariants
        envs = {
            "chrom_col": "chr",
            "pos_col": "pos",
            "pval_col": args.var_pval_col,
            # "label_col": "SNP",
            "chroms": args.chroms,
            "zoom": args.chroms,
            "signif": args.signif,
            # "title": f"Manhattan Plot (Combined {args.pval_col})",
            "ylabel": f"-log10(Combined {args.var_pval_col})",
        }

    class ManhattanPlots(Pipen):
        if args.ceqtl_vars:
            starts = AddCoordsToTrios, AddCoordsToVariants
        else:
            starts = AddCoordsToTrios

    return ManhattanPlots()


def main():
    args = parser.parse_extra_args()
    try:
        pipen = get_pipeline(args)
    except AttributeError as e:
        print("Error in arguments: {}".format(e))
        parser.print_help()
    else:
        pipen.run()
