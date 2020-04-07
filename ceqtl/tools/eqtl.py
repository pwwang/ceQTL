"""cis-eQTL calling"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.eqtl import pMatrixeQTL
from bioprocs.tsv import pTsvHeader, pTsvJoin, pTsvColSelect
from bioprocs.common import pSort

def pipeline(opts):
    """Construct the pipeline"""
    # we need get the intersect samples for both expr and gtype
    pTsvHeader.input = [opts.gtype, opts.expr]

    pSort.depends = pTsvHeader
    pSort.args.inopts.skip = 1

    pTsvJoin.depends = pSort
    pTsvJoin.input = lambda ch: [ch.flatten()]
    pTsvJoin.args.inopts.cnames = True
    # also add ROWNAME, opts.expr does not have it
    pTsvJoin.args.helper = [
        'rowname_added = False',
        'def write(writer, r1, r2):',
        '   global rowname_added',
        '   if not rowname_added:',
        '       writer.write(["ROWNAME"])',
        '       writer.write(["ID"])',
        '       rowname_added = True',
        '   writer.write(r1)'
    ]
    pTsvJoin.args.do = 'lambda writer, r1, r2: write(writer, r1, r2)'

    pTsvColSelect.depends = pTsvJoin
    pTsvColSelect.input = lambda ch: ch.rep_row(2).insert(0, [opts.gtype, opts.expr])

    pMatrixeQTL.depends = pTsvColSelect
    pMatrixeQTL.input = lambda ch: ch.unfold()
    pMatrixeQTL.args.pval = opts.pval
    pMatrixeQTL.args.snppos = opts.snppos
    pMatrixeQTL.args.genepos = opts.genepos
    pMatrixeQTL.args.dist = opts.dist
    pMatrixeQTL.args.model = 'modelANOVA'
    pMatrixeQTL.output = ('outfile:file:%s,'
                          'alleqtl:file:{{i.snpfile | fn}}-{{i.expfile | fn}}.eqtl.txt'
                          ) % Path(opts.outfile).name
    pMatrixeQTL.config.export_dir = Path(opts.outfile).parent
    pMatrixeQTL.config.export_part = "outfile"

    return pTsvHeader

def main(opts):
    """Main function"""

    PyPPL(name='cis-eQTL calling') \
        .start(pipeline(opts)) \
        .run()
