"""Implementation of ceQTL-tools addtfgenes"""
from pathlib import Path
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import logger

def read_tfgenes(tgfile):
    """Read tf-gene pairs"""
    logger.info('Reading TF-gene pairs in %s ...', tgfile)
    reader = TsvReader(tgfile, cnames=False)
    ret = {} # gene => tf
    for row in reader:
        ret.setdefault(row[1], set()).add(row[0])
    reader.close()
    return ret

def main(opts):
    """Main function"""
    org_tfgenes = read_tfgenes(opts.origin)
    add_tfgenes = read_tfgenes(opts.addition)
    writer = TsvWriter(opts.outfile)
    logger.info('Writing the union set to %s ...', opts.outfile)
    for gene, tfs in org_tfgenes.items():
        for tf in (tfs | add_tfgenes.pop(gene, set())):
            writer.write([tf, gene])
    for gene, tfs in add_tfgenes.items():
        for tf in tfs:
            writer.write([tf, gene])
    writer.close()
    logger.info('Done.')
