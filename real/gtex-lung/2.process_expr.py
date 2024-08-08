"""Process the expression from GTex"""
# biopipen v0.29.0
from pipen import Pipen
from pipen_args import parser
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.misc import Shell
from biopipen.ns.gene import GeneNameConversion

parser.add_argument("--tfgmt", help="TF-gene GMT file", required=True)
parser.add_argument("--expr", help="Expression file", required=True)
args = parser.parse_extra_args()


class ConvertToMatrix(Shell):
    """Convert the expression to a matrix"""
    input_data = [args.expr]
    output = "outfile:file:{{in.infile | stem | stem}}.matrix.txt"
    envs = {
        "cmd": "zcat $infile | cut -d $'\\t' -f 4,5- > $outfile"
    }


class ConvertToGeneSymbol(GeneNameConversion):
    requires = ConvertToMatrix
    envs = {"infmt": "ensembl.gene", "notfound": "skip", "output": "replace"}


# Collect all genes from TF-gene pairs
class ExtractAllGenes(Shell):
    input_data = [args.tfgmt]
    envs = {"cmd": "cut -f1,3- $infile | tr '\t' '\n' | sort -u > $outfile"}


class FilterGenes(Proc):
    requires = ExtractAllGenes, ConvertToGeneSymbol
    input = "agfile:file, mtfile:file"
    output = "outfile:file:{{in.mtfile | basename}}"
    lang = config.lang.rscript
    script = """
        agfile <- {{in.agfile | r}}
        mtfile <- {{in.mtfile | r}}
        outfile <- {{out.outfile | r}}

        mat <- read.table(
            mtfile,
            header = TRUE,
            row.names = NULL,
            sep = "\\t",
            check.names = FALSE
        )
        genes <- readLines(agfile)
        matgenes <- mat[, 1, drop = TRUE]
        genes <- genes[genes %in% matgenes]
        mat <- mat[matgenes %in% genes, , drop = FALSE]
        genes <- mat[, 1, drop = TRUE]
        mat <- mat[!duplicated(genes), , drop = FALSE]
        genes <- mat[, 1, drop = TRUE]
        mat <- mat[, -1, drop = FALSE]
        rownames(mat) <- genes
        write.table(
            mat,
            outfile,
            col.names = TRUE,
            row.names = TRUE,
            sep = "\\t",
            quote = FALSE
        )
    """


class ProcessExpr(Pipen):
    __doc__ = __doc__
    starts = ConvertToMatrix, ExtractAllGenes


if __name__ == "__main__":
    ProcessExpr().run()
