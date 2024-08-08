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
        in_score: The column in the input file that are scores
        in_dir: If `-`, the score will be negated
        gold_ids (list): The columns in the gold file that are ids
        gold_score: The column in the gold file that are scores
        gold_cutoff (type=float): The cutoff on gold_score to be positive
        gold_dir (choice): The direction of gold_score to be positive
            +: larger is positive
            -: smaller is positive
        gold_uniq: Make sure the gold ids are unique
        gold_score_agg: Aggregate the gold scores if ids are not unique.
            An R function that takes the scores as input and return a single score.
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
        "gold_uniq": False,
        "gold_score_agg": "sumlog",
    }
    script = """
        library(rlang)
        library(dplyr)
        library(tidyr)
        library(metap)

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
        gold_uniq <- {{envs.gold_uniq | r}}
        gold_score_agg <- {{envs.gold_score_agg}}

        agg <- function(scores) {
            if (length(scores) == 1) {
                return(scores)
            }
            return(gold_score_agg(scores))
        }

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
            if (gold_uniq) {
                gold <- gold %>%
                    group_by(!!!syms(gold_ids)) %>%
                    summarise(!!sym(gold_score) := agg(!!sym(gold_score)))
            }
            if (!is.null(gold_cutoff)) {
                if (gold_dir == "+") {
                    # gold <- gold[gold[[gold_score]] >= gold_cutoff, ]
                    gold$gold_label <- gold[[gold_score]] >= gold_cutoff
                } else {
                    # gold <- gold[gold[[gold_score]] <= gold_cutoff, ]
                    gold$gold_label <- gold[[gold_score]] <= gold_cutoff
                }
            } else {
                gold$gold_label <- TRUE
            }
            gold$label <- gold$gold_label
            gold$gold_label <- NULL
            gold[[gold_score]] <- NULL
        }
        # get all the positive ids
        # golds <- gold[, gold_ids]
        # concat all rows
        # golds <- apply(golds, 1, paste, collapse = ".")
        gold <- unite(gold, "id", all_of(gold_ids), sep = ".") %>% select(id, label)

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
        # in_labels <- in_records %in% golds
        in_scores <- indata[[in_score]]
        if (in_dir == "-") {
            in_scores <- -in_scores
        }
        # out <- data.frame(id = in_records, in_label = in_labels, in_score = in_scores)
        out <- data.frame(id = in_records, score = in_scores)
        out <- full_join(out, gold, by = "id")
        minscore <- min(out$score, na.rm = T)
        out$score[is.na(out$score)] <- minscore - runif(sum(is.na(out$score))) * 1e-6
        out$label[is.na(out$label)] <- FALSE
        out <- out %>% select(id, label, score)
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
