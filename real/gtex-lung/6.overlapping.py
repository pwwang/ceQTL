"""Generate venn diagram for overlapping variants/variant-gene pairs between
eqtls and ceqtls."""

from pipen import Pipen
from pipen_args import parser
from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.plot import VennDiagram as VennDiagramPlot


parser.add_argument(
    '--eqtls',
    help='eqtl file',
    required=True,
)
parser.add_argument(
    '--ceqtls',
    help='ceqtl file (variant-gene pairs)',
    required=True,
)
parser.add_argument(
    '--eqtl-cols',
    help='columns to use for eqtls identity',
    default="SNP,gene",
)
parser.add_argument(
    '--ceqtl-cols',
    help='columns to use for ceqtls identity',
    default="SNP,Target",
)
parser.add_argument(
    '--eqtl-pos',
    help='expression to select positive eqtls',
    default="FDR < 1e-10",
)
parser.add_argument(
    '--ceqtl-pos',
    help='expression to select positive ceqtls',
    default="MetaPadj < 1e-10",
)


def pipeline(args):
    class PrepareVennDiagram(Proc):
        input = "efile:file,cefile:file"
        input_data = [(args.eqtls, args.ceqtls)]
        output = "outfile:file:venn_diagram.txt"
        lang = config.lang.rscript
        envs = {
            "eqtl_cols": args.eqtl_cols,
            "ceqtl_cols": args.ceqtl_cols,
            "eqtl_pos": args.eqtl_pos,
            "ceqtl_pos": args.ceqtl_pos,
        }
        export = True
        script = """
            library(dplyr)
            library(tidyr)

            efile = {{in.efile | r}}
            cefile = {{in.cefile | r}}
            outfile = {{out.outfile | r}}
            eqtl_cols = {{envs.eqtl_cols | r}}
            ceqtl_cols = {{envs.ceqtl_cols | r}}
            eqtl_cols = strsplit(eqtl_cols, ",")[[1]]
            ceqtl_cols = strsplit(ceqtl_cols, ",")[[1]]

            eqtls = read.table(efile, header=TRUE, sep="\t")
            ceqtls = read.table(cefile, header=TRUE, sep="\t")
            eqtls = eqtls %>% filter({{envs.eqtl_pos}}) %>%
                select(all_of(eqtl_cols)) %>%
                unite("eQTL", everything(), sep="-") %>%
                pull(eQTL) %>%
                unique()
            ceqtls = ceqtls %>% filter({{envs.ceqtl_pos}}) %>%
                select(all_of(ceqtl_cols)) %>%
                unite("ceQTL", everything(), sep="-") %>%
                pull(ceQTL) %>%
                unique()
            all = union(eqtls, ceqtls) %>% unique()
            df = data.frame(
                eQTL = all %in% eqtls,
                ceQTL = all %in% ceqtls
            )
            rownames(df) = all
            write.table(df, outfile, sep="\\t", quote=FALSE)
        """

    class VennDiagram(VennDiagramPlot):
        requires = PrepareVennDiagram
        envs = {
            "intype": "computed",
            "inopts": {"row.names": 1, "header": True},
            "args": {"edge_size": 0},
            "ggs": ["scale_fill_distiller(palette = 'YlGnBu', direction = 1)"],
        }

    class eQTLceQTLOverlapping(Pipen):
        starts = PrepareVennDiagram

    return eQTLceQTLOverlapping()


if __name__ == "__main__":
    args = parser.parse_extra_args()
    pipeline(args).run()
