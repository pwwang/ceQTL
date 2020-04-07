"""Custom processes using by ceQTL pipeline"""

from diot import Diot
from bioprocs import proc_factory, params

pTFT2GeneGroups = proc_factory(
    desc='Convert TF-target relations to row group file',
    input='infile:file',
    output='outfile:file:{{i.infile | stem2}}.rg.txt',
    lang=params.Rscript.value,
    script='file:scripts/pTFT2GeneGroups.R'
)

pSGPair2Beds = proc_factory(
    desc='Convert SNP-gene pairs to BED files for simulation',
    lang=params.python.value,
    input='infile:file',
    output=('snpfile:file:{{i.infile | stem2}}.snps.bed, '
            'genefile:file:{{i.infile | stem2}}.gene.bed'),
    args=Diot(snppergene=5, nchr=10, seed=0, dist=1000),
    script='file:scripts/pSGPair2Beds.py'
)

pTFTSG2MedCases = proc_factory(
    desc='Convert TF-target and snp-gene pair file to mediation case file',
    lang=params.python.value,
    input='tftfile:file, sgfile:file',
    output='outfile:file:{{i.tftfile | fn}}.medcase.txt',
    script='file:scripts/pTFTSG2MedCases.py'
)
