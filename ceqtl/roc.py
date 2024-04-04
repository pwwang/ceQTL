from pipen import Pipen
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import ROC as _ROC


class Input(Proc):
    """Prepare input for ROC plot

    Input:
        infile: The input file
        goldfile: The gold standard file

    Output:
        outfile: The output file

    Envs:
        in_ids (list): The columns in the input file that are ids
        gold_ids (list): The columns in the gold file that are ids
        in_score: The column in the input file that are scores
        gold_score: The column in the gold file that are scores
        gold_cutoff (type=float): The cutoff on gold_score to be positive
        gold_dir (choice): The direction of gold_score to be positive
            +: larger is positive
            -: smaller is positive
    """
    input = "infile:file, goldfile:file"
    output = "outfile:file:{{in.infile | stem}}.4roc.txt"
    lang = config.lang.rscript
    envs = {
        "in_ids": [1],
        "in_score": [2],
        "in_dir": "+",  # larger is positive
        "gold_ids": [1],
        "gold_score": None,  # all are positive
        "gold_cutoff": None,  # The cutoff on gold_score to be positive
        "gold_dir": "+",  # larger is positive
    }
    script = """
        infile <- {{in.infile | r}}
        goldfile <- {{in.goldfile | r}}
        outfile <- {{out.outfile | r}}
        in_ids <- {{envs.in_ids | r}}
        in_score <- {{envs.in_score | r}}
        in_dir <- {{envs.in_dir | r}}
        gold_ids <- {{envs.gold_ids | r}}
        gold_score <- {{envs.gold_score | r}}
        gold_cutoff <- {{envs.gold_cutoff | r}}
        gold_dir <- {{envs.gold_dir | r}}

        # read the gold standard
        gold <- read.table(
            goldfile, header = T, stringsAsFactors = F,
            row.names = NULL, check.names = F, sep = "\t")
        if (is.numeric(gold_ids)) {
            gold_ids <- colnames(gold)[gold_ids]
        }
        if (is.numeric(gold_score)) {
            gold_score <- colnames(gold)[gold_score]
        }
        if (!is.null(gold_score)) {
            if (gold_dir == "+") {
                gold <- gold[gold[[gold_score]] >= gold_cutoff, ]
            } else {
                gold <- gold[gold[[gold_score]] <= gold_cutoff, ]
            }
        }
        # get all the positive ids
        golds <- gold[, gold_ids]
        # concat all rows
        golds <- apply(golds, 1, paste, collapse = ".")

        # read the input
        indata <- read.table(
            infile, header = T, stringsAsFactors = F,
            row.names = NULL, check.names = F, sep = "\t")
        if (is.numeric(in_ids)) {
            in_ids <- colnames(indata)[in_ids]
        }
        if (is.numeric(in_score)) {
            in_score <- colnames(indata)[in_score]
        }
        in_records <- apply(indata[in_ids], 1, paste, collapse = ".")
        in_labels <- in_records %in% golds
        in_scores <- indata[[in_score]]
        if (in_dir == "-") {
            in_scores <- -in_scores
        }
        out <- data.frame(id = in_records, label = in_labels, score = in_scores)
        write.table(out, file = outfile, quote = F, sep = "\t", row.names = F)
    """


class ROC(_ROC):
    """"""
    requires = Input
    envs = {"fdr": True}


class Pipeline(Pipen):
    name = "ceqtl_roc"
    desc = "Running ROC analysis"
    starts = Input


def main():
    Pipeline().run()
