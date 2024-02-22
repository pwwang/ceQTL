
from biopipen.core.proc import Proc
from biopipen.core.config import config


class DataPreparation(Proc):
    """Prepare the data for the analysis.

    This process:
    1. checkes whether the number and order of samples in the genotype matrix and the
    expression matrix are the same.
    2. checks whether the number and order of samples in the covariate matrix and the
    expression matrix are the same.
    2. splits the SNPs into chunks to run in parallel.
    3. saves the split expression + covariate matrix, genotype matrix and
    the formula file for each chunk.

    Input:
        geno: The genotype matrix file, with SNP as columns and samples as rows.
        expr: The expression matrix file
        cov: The covariate variables matrix.
        snpgene: The SNP-gene file in GMT format
        tftarget: The TF-target file in GMT format

    Output:
        outdir: The output directory

    Envs:
        nchunks: The number of chunks to break the SNPs into.
        transpose_geno (flag): Whether to transpose the genotype matrix if
            the SNPs are in rows.
        transpose_expr (flag): Whether to transpose the expression matrix if
            the genes are in rows.
        transpose_cov (flag): Whether to transpose the covariate matrix if
            the covariates are in rows.
    """
    input = "geno:file, expr:file, cov:file, snpgene:file, tftarget:file"
    output = "outdir:dir:{{in.geno | stem}}.prepared"
    lang = config.lang.rscript
    envs = {
        "nchunks": 4,
        "transpose_geno": False,
        "transpose_expr": False,
        "transpose_cov": False,
    }
    script = "file://scripts/DataPreparation.R"
