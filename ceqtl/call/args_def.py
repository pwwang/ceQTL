"""Provide function to add extra arguments to the parser"""
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pipen_args import Parser


def add_args(parser: Parser) -> Parser:
    """Add extra arguments to the parser"""
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
        "--tftarget",
        help="The TF-target file in GMT format",
        required=True,
        type="path",
    )
    parser.add_argument(
        "--genesnp",
        help="The Gene-SNP file in GMT format",
        required=True,
        type="path",
    )
    parser.add_argument(
        "--nchunks",
        help="Break the SNPs into chunks to run in parallel",
        required=True,
        type="path",
    )
    parser.add_argument(
        "--ncores",
        help="Number of cores to use for all processes",
        default=1,
    )
    return parser
