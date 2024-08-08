from pipen import Pipen
from biopipen.ns.snp import MatrixEQTL as _MatrixEQTL


class MatrixEQTL(_MatrixEQTL):
    output = "alleqtls:file:alleqtls.txt, cisqtls:file:cisqtls.txt"
    envs = {"fdr": True, "pval": 1}


class Pipeline(Pipen):
    name = "ceqtl_eqtl"
    desc = "Running eQTL analysis"
    starts = MatrixEQTL
    plugin_opts = {"args_flatten": False}


def main():
    Pipeline().run()
