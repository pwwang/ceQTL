"""Plot ROC curves for different measurements"""

# pylint: disable=invalid-name,assigning-non-slot
from pathlib import Path
from pyppl import PyPPL
from bioprocs.tsv import pTsvColSelect, pTsvJoin, pTsvHeader, pTsv, pTsvCbind, pTsvReplaceHeader
from bioprocs.common import pSort, pShell
from bioprocs.plot import pROC

def colselect(infile, col, cut, rev, tag):
    """Select the columns from input file"""
    pTsvColSelect1 = pTsvColSelect.copy(tag=tag)
    if str(col).isdigit():
        pTsvColSelect1.input = [infile]
        pTsvColSelect1.args.cols = [0, int(col)]
        start = end = pTsvColSelect1
    else:
        pTsvHeader1 = pTsvHeader.copy(tag=tag)
        pTsvHeader1.input = [infile]
        start = pTsvHeader1

        pShell1 = pShell.copy(tag=tag)
        pShell1.depends = pTsvHeader1
        pShell1.args.cmd = """
        head -1 {{i.infile | quote}} > {{o.outfile | quote}}
        echo "%s" >> {{o.outfile | quote}}
        """ % col

        pTsvColSelect1.depends = pShell1
        pTsvColSelect1.input = 'infile:file, colfile:file'
        pTsvColSelect1.input = lambda ch: ch.insert(0, infile)
        start = pTsvHeader1
        end = pTsvColSelect1

    if cut and cut < 1:
        pTsv1 = pTsv.copy(tag=tag)
        pTsv1.depends = pTsvColSelect1
        pTsv1.args.row = 'lambda row: float(row[%r]) < %f' % (col, cut)
        end = pTsv1
    elif cut and cut > 1:
        pSort1 = pSort.copy(tag=tag)
        pSort1.depends = pTsvColSelect1
        pSort1.args.inopts.skip = 1
        pSort1.args.params.k = '2g%s' % ('' if rev else 'r')

        pShell2 = pShell.copy(tag=tag)
        pShell2.depends = pSort1
        pShell2.args.cmd = "head -n %d {{i.infile | quote}} > {{o.outfile | quote}}" % (int(cut)+1)
        end = pShell2

    if rev:
        pTsv2 = pTsv.copy(tag=tag)
        pTsv2.depends = end
        pTsv2.args.row = 'lambda row: row.__setitem__({0!r}, "-" + row[{0!r}])'.format(col)
        end = pTsv2
    return start, end

def pipeline(opts):
    """Construct the pipeline"""

    starts = []
    ends = []
    for i, infile in enumerate(opts.infiles):
        start, end = colselect(infile, opts.cols[i],
                               opts.cuts[i] if opts.cuts else None,
                               opts.rev, "in%s" % (i+1))
        starts.append(start)
        ends.append(end)

    # duplicate gold standard columns
    pTsvGold = pTsv.copy()
    pTsvGold.input = [opts.gold]
    pTsvGold.args.inopts.cnames = False
    pTsvGold.args.row = 'lambda row: [row[0], 1]'
    starts.append(pTsvGold)

    pROC.args.params.bestCut = False
    pROC.args.ggs.theme_bw = {}

    if not opts.sep:
        pSort.depends = ends
        pSort.input = lambda *chs: [ch.get() for ch in chs]
        pSort.args.inopts.skip = 1

        pTsvJoin.depends = pSort
        pTsvJoin.input = lambda ch: [ch.flatten()]
        #pTsvJoin.args.outopts.cnames = ["ROWNAME"] + [Path(infile).stem for infile in opts.infiles]
        pTsvJoin.args.inopts.cnames = True
        pTsvJoin.args.do = "lambda writer, *rs: writer.write([rs[0][0]] + [r[1] for r in rs])"

        pTsvCbind.depends = pTsvGold, pTsvJoin
        pTsvCbind.input = lambda ch1, ch2: [ch1.cbind(ch2).flatten()]
        pTsvCbind.args.inopts.cnames = False
        pTsvCbind.args.inopts.rnames = True
        pTsvCbind.args.fill = True
        pTsvCbind.args.fn2cname = "NULL"

        # set the non-hit record to 0
        # and remove NAs from predictions
        pTsvGoldFalse = pTsv.copy()
        pTsvGoldFalse.depends = pTsvCbind
        pTsvGoldFalse.args.row = ('lambda row: False '
                                  'if any(r == "NA" for r in row[2:]) '
                                  'else row.__setitem__(1, 0) '
                                  'if row[1] == "NA" '
                                  'else None')

        pTsvReplaceHeader.depends = pTsvGoldFalse
        pTsvReplaceHeader.args.cnames = ["ROWNAME", "GOLD"] + opts.names

        pROC.depends = pTsvReplaceHeader
    else:
        # add header
        pTsvReplaceHeader.depends = pTsvGold
        pTsvReplaceHeader.args.inopts.cnames = False
        pTsvReplaceHeader.args.cnames = ['ROWNAME', 'GOLD']

        ends.insert(0, pTsvReplaceHeader)
        pTsvCbind.depends = ends
        pTsvCbind.output = 'outfile:file:{{i.infiles | [-1] | stem2 }}.cbound.txt'
        pTsvCbind.input = (lambda ch_gold, *chs: [ch_gold.cbind(ch).flatten() for ch in chs])

        pTsvGoldFalse = pTsv.copy()
        pTsvGoldFalse.depends = pTsvCbind
        pTsvGoldFalse.config.export_dir = opts.outdir
        pTsvGoldFalse.args.row = ('lambda row: False '
                                  'if any(r == "NA" for r in row[2:]) '
                                  'else row.__setitem__(1, 0) '
                                  'if row[1] == "NA" '
                                  'else None')

        pROC.depends = pTsvGoldFalse
        pROC.config.export_dir = opts.outdir

    return starts

def main(opts):
    """Main function"""

    PyPPL(forks=opts.nthread).start(pipeline(opts)).run()
