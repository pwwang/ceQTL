"""Implementation of ceQTL-tools genotype-qc"""
# pylint: disable=invalid-name,assigning-non-slot
from pyppl import PyPPL
from bioprocs.tsv import pMatrixR
from bioprocs.rnaseq import pUnitConversion

def median_to_rfunc(median):
    return '''
    function(expr) {
        meds = apply(expr, 1, median)
        expr[meds > %s, ]
    }
    ''' % median

def pipeline(opts):
    """Construct the pipeline"""
    start_processes = []
    end_process = None

    if opts.premed or opts.inform:
        pPreMedianFilter = pMatrixR.copy()
        pPreMedianFilter.input = [opts.expr]
        pPreMedianFilter.args.inopts.dup = 'mean'
        pPreMedianFilter.args.code = []
        if opts.premed:
            pPreMedianFilter.args.code.append(
                "mat = (%s)(mat)" % median_to_rfunc(opts.premed).strip()
            )
        if opts.inform:
            pPreMedianFilter.args.code.append(
                "mat = (%s)(mat)" % opts.inform.strip()
            )

        start_processes.append(pPreMedianFilter)
        end_process = pPreMedianFilter

    if opts.outunit:
        if start_processes:
            pUnitConversion.depends = start_processes
        else:
            pUnitConversion.input = [opts.expr]
            start_processes.append(pUnitConversion)
        pUnitConversion.args.inopts.dup = 'mean'
        pUnitConversion.args.inunit = opts.inunit
        pUnitConversion.args.outunit = opts.outunit
        end_process = pUnitConversion

    if opts.postmed or opts.outform:
        pPostMedianFilter = pMatrixR.copy()
        pPostMedianFilter.args.inopts.dup = 'mean'
        pPostMedianFilter.args.code = []
        if start_processes:
            pPostMedianFilter.depends = start_processes
        else:
            pPostMedianFilter.input = [opts.expr]
            start_processes.append(pPostMedianFilter)

        if opts.outform:
            pPostMedianFilter.args.code.append(
                "mat = (%s)(mat)" % opts.outform.strip()
            )
        if opts.postmed:
            pPostMedianFilter.args.code.append(
                "mat = (%s)(mat)" % median_to_rfunc(opts.postmed).strip()
            )
        end_process = pPostMedianFilter

    end_process.config.export_dir = opts.outdir

    return start_processes

def main(opts):
    """Main function"""

    if opts.outform == 'log2':
        opts.outform = 'function(x) log2(x+1)'

    PyPPL(name='Quality control for expression data') \
        .start(pipeline(opts)) \
        .run()
