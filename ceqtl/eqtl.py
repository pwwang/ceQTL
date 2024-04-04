from pipen import Pipen
from biopipen.ns.snp import MatrixEQTL as _MatrixEQTL


class MatrixEQTL(_MatrixEQTL):
    envs = {"fdr": True}


class Pipeline(Pipen):
    name = "ceqtl_eqtl"
    desc = "Running eQTL analysis"
    starts = MatrixEQTL


def main():
    Pipeline().run()
