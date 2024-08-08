source("{{biopipen_dir}}/utils/misc.R")

# library(foreach)
# library(doParallel)
# library(nleqslv)
library(ggplot2)

genofile <- {{in.genofile | r}}
seed <- {{in.seed | r}}
exfile <- {{out.exfile | r}}
snpgenefile <- {{out.snpgenefile | r}}
tftargetfile <- {{out.tftargetfile | r}}
truthallfile <- {{out.truthallfile | r}}
truthfile <- {{out.truthfile | r}}
truthplotdir <- {{out.truthfile | stem | append: ".plots" | r}}
dir.create(truthplotdir, showWarnings = FALSE)

ncores <- {{envs.ncores | r}}
ngenes <- {{envs.ngenes | r}}
npair_range <- {{envs.npair_range | r}}
nregulator_range <- {{envs.nregulator_range | r}}
ntarget_range <- {{envs.ntarget_range | r}}
gene_prefix <- {{envs.gene_prefix | r}}
gene_index_start <- {{envs.gene_index_start | r}}
noise <- {{envs.noise | r}}
transpose_exprs <- {{envs.transpose_exprs | r}}

set.seed(seed)
# registerDoParallel(cores=24)

#### Generate genotype matrix
suppressPackageStartupMessages({
    library(Rfast)
    library(optimg)
    library(parallel)
})

#' Generate a group matrix, with rows representing cases and
#' columns representing samples.
#'
#' @param x An integer, number of cases to generate
#' @param nsamples An integer, number of samples to generate
#' @param ngroups_range An integer vector of length 1 or 2, range of number
#'   of groups per each case
#' @param group_minsample An integer, minimum number of samples per group
#' @param sample_prefix A string, prefix for sample names
#' @param sample_index_start An integer, starting index for sample names
#' @param case_prefix A string, prefix for case names
#' @param case_index_prefix An integer, starting index for case names
#' @return A data frame, with rows representing cases and columns representing
#' @export
as_groupmat.int <- function(
    x,
    nsamples,
    ngroups_range,
    group_minsample = 5,
    sample_prefix = "Sample",
    sample_index_start = 1,
    case_prefix = "Case",
    case_index_prefix = 1,
    ...
) {
    cases <- list()
    if (length(ngroups_range) == 1) { ngroups_range <- rep(ngroups_range, 2) }
    minprob <- group_minsample / nsamples
    for (i in seq_len(x)) {
        casename <- paste0(case_prefix, case_index_prefix + i - 1)
        if (ngroups_range[1] == ngroups_range[2]) {
            ngroups <- ngroups_range[1]
        } else {
            ngroups <- sample(x = ngroups_range[1]:ngroups_range[2], size = 1)
        }
        if (ngroups * group_minsample > nsamples) {
            stop("ngroups * group_minsample > nsamples")
        }
        probs <- 0
        while (min(probs) < minprob) {
            probs <- rnorm(ngroups)
            probs <- probs / sum(probs)
        }
        case <- list(ngroups = ngroups, probs = probs)
        cases[[casename]] <- case
    }
    as_groupmat.list(
        cases,
        nsamples,
        sample_prefix = sample_prefix,
        sample_index_start = sample_index_start
    )
}

#' @param x A list of cases, each case is a list with elements ngroups and
#'   probs to draw samples for each group
#' @param nsamples An integer, number of samples to generate
#' @param probs_fn A function or string, function to generate probabilities
#'   for each group, or a string to be evaluated as a function
#' @param group_minsample An integer, minimum number of samples per group
#' @param sample_prefix A string, prefix for sample names
#' @param sample_index_start An integer, starting index for sample names
#' @return A data frame, with rows representing cases and columns representing
#' @export
as_groupmat.list <- function(
    x,
    nsamples,
    probs_fn = rnorm,
    group_minsample = 5,
    sample_prefix = "Sample",
    sample_index_start = 1,
    ...
) {
    if (is.character(probs_fn)) { probs_fn <- eval(parse(text = probs_fn)) }
    minprob <- group_minsample / nsamples
    samples <- paste0(sample_prefix, (1:nsamples) + sample_index_start - 1)
    groupmat <- c()
    for (n in names(x)) {
        case <- x[[n]]
        ngroups <- case$ngroups
        probs <- case$probs
        if (is.null(probs)) {
            probs <- 0
            while (min(probs) < minprob) {
                probs <- probs_fn(ngroups)
                probs <- probs / sum(probs)
            }
        } else if (length(probs) == ngroups - 1) {
            probs <- c(probs, 1 - sum(probs))
        }
        groups <- sample(1:ngroups, size = nsamples, prob = probs, replace = TRUE)
        groupmat <- rbind(groupmat, groups)
        rownames(groupmat)[nrow(groupmat)] = n
    }
    groupmat <- as.data.frame(groupmat)
    colnames(groupmat) <- samples

    return(groupmat)
}

#' @param x A data frame
#' @return The input data frame
#' @export
as_groupmat.data.frame <- function(x, ...) { x }

invisible(setGeneric("as_groupmat", function(x, ...) standardGeneric("as_groupmat")))
setMethod("as_groupmat", "list", as_groupmat.list)
setMethod("as_groupmat", "integer", as_groupmat.int)
setMethod("as_groupmat", "numeric", as_groupmat.int)
setMethod("as_groupmat", "data.frame", as_groupmat.data.frame)
### End of as_groupmat


#' Generate expected correlations for each case and gene pairs
#'
#' @param groupmat A group matrix. Columns are samples and rows are cases.
#'   The values are usually starting from 1 and represent the group of the samples.
#' @param npair_range A vector of length 2. The range of the number of gene pairs
#'   for each case.
#' @param gene_pairs A data frame with rows as regulators and columns as targets.
#'   The values are binary indicating whether the regulator regulates the target.
#' @param mincordiff A numeric, the minimum difference between the expected
#' @return A list of lists. The outer list has the same names as the rows of groupmat,
#'   which is the case names. The inner list is a list of data frames. The names of
#'   the inner list are sample groups. The data frames have the same dimension as
#'   gene_pairs. The values are the expected correlations between the regulators and
#'   the targets.
#' @export
generate_params <- function(groupmat, npair_range, gene_pairs, mincordiff = 0.5) {
    find_unmet_idx <- function(corvecs, idxes) {
        combined <- combn(names(corvecs), 2)
        diff <- NULL
        for (i in seq_len(ncol(combined))) {
            pair <- combined[, i]
            if (is.null(diff)) {
                diff <- abs(corvecs[[pair[1]]] - corvecs[[pair[2]]])
            } else {
                # element-wise max
                diff <- pmax(diff, abs(corvecs[[pair[1]]] - corvecs[[pair[2]]]))
            }
        }
        return (seq_along(diff)[diff < mincordiff])
    }
    out = list()
    # list(
    #    c1 = list(
    #       1 = cordf1,
    #       2 = cordf2,
    #       3 = cordf3
    #    ),
    #    c2 = ...
    # )
    #
    # cordf:
    #      g1  g2  g3
    # tf1  .1  .2  .3
    # tf2  .4  .5  .6
    # ...
    allgps <- as.vector(as.matrix(gene_pairs))
    # The edges of the network
    allindexes <- seq_along(allgps)[allgps == 1]
    for (cs in rownames(groupmat)) {
        ugroups <- as.character(sort(unique(unlist(groupmat[cs, ]))))
        npairs <- sample(npair_range[1]:npair_range[2], 1)

        # pick the edges for the case
        pair_idxes <- sample(allindexes, npairs)

        # generate random correlations
        out[[cs]] = lapply(ugroups, function(g) runif(length(pair_idxes), -1 , 1))
        names(out[[cs]]) <- ugroups
        idxes <- find_unmet_idx(out[[cs]], seq_len(pair_idxes))
        while (length(idxes) > 0) {
            for (g in names(out[[cs]])) {
                out[[cs]][[g]][idxes] <- runif(length(idxes), -1, 1)
            }
            idxes <- find_unmet_idx(out[[cs]], idxes)
        }
        out[[cs]] <- lapply(out[[cs]], function(x) {
            corvecs <- rep(NA, length(allgps))
            corvecs[pair_idxes] <- x
            corvecs <- matrix(corvecs, nrow = nrow(gene_pairs))
            colnames(corvecs) <- colnames(gene_pairs)
            rownames(corvecs) <- rownames(gene_pairs)
            # remove all columns with NAs
            corvecs <- corvecs[, colSums(is.na(corvecs)) < nrow(corvecs), drop = FALSE]
            # remove all rows with NAs
            corvecs <- corvecs[rowSums(is.na(corvecs)) < ncol(corvecs), , drop = FALSE]
            return(corvecs)
        })
    }

    return(out)
}

#' Convert gene pairs data frame to a data frame that can be saved as a GMT file
#'
#' @param gene_pairs A dataframe with rows as regulators and columns as targets.
#'   The values are binary indicating whether the regulator regulates the target.
#' @return A data frame with 3 columns. The first column is the regulator, the second
#'   column is empty, and the third column is a string of targets separated by tabs.
#' @export
genepairs_to_gmt <- function(gene_pairs) {
    # convert to gmt dataframe
    # gene_pairs = c("g1 -> g2", "g1 -> g3", "g2 -> g3")
    out <- c()
    for (reg in rownames(gene_pairs)) {
        targets <- colnames(gene_pairs[, gene_pairs[reg, ] == 1, drop = FALSE])
        out <- rbind(out, data.frame(
            V1 = reg,
            V2 = "",
            V3 = paste(targets, collapse = "\t")
        ))
    }
    return(out)
}

#' Convert params to a data frame that can be saved as a GMT file
#'
#' @param params A list of lists. The outer list has the same names as the rows of
#'   groupmat, which is the case names. The inner list is a list of data frames. The
#'   names of the inner list are sample groups.
#' @return A data frame with 3 columns. The first column is the case name, the second
#'   column is empty, and the third column is a string of regulators separated by tabs.
params_to_gmt <- function(params) {
    # convert to gmt dataframe
    # params = list(
    #    c1 = list(`1` => cordf1, `2` => cordf2, `3` => cordf3)
    # )
    # cordf:
    #     g1  g2  g3
    # tf1 NA  .2  .3
    # tf2 .4  .5  NA
    out <- c()
    for (cs in names(params)) {
        targets <- colnames(params[[cs]][["1"]])
        out <- rbind(out, data.frame(
            V1 = cs,
            V2 = "",
            V3 = paste(targets, collapse = "\t")
        ))
    }
    return (out)
}


#' The cost function to minimize
#'
#' @param weights A numeric vector, the weights for the expression matrix
#' @param x A numeric matrix, the expression matrix
#' @param corrs A list of lists. The outer list has the same names as the rows of
#'   groupmat, which is the case names. The inner list is a list of data frames. The
#'   names of the inner list are sample groups.
#' @param groupmat A group matrix. Columns are samples and rows are cases.
#'  The values are usually starting from 1 and represent the group of the samples.
#' @param costfn A string, the cost function to use. It can be "ew" for element-wise
#'   cost, or "cor" for correlation cost.
#' @param ncores An integer, the number of cores to use for parallel computing
#' @return A numeric, the cost
cost_fn <- function(weights, x, corrs, groupmat, costfn, ncores){
    # corrs = list(
    #    c1 = list( 1 => cordf1, 2 => cordf2, 3 => cordf3)
    # )
    # cordf:
    #     g1  g2  g3
    # tf1 NA  .2  .3
    # tf2 .4  .5  NA
    # x:
    #    g1 g2 g3
    # s1 .1 .2 .3
    # s2 .4 .5 .6
    # groupmat:
    #       s1	s2	s3
    # Case1	3	1	1
    # Case2	1	3	2
    # Case3	1	3	3
    x <- weights * x
    # Spread the case and its groups for parallel computing
    args <- list()
    for (cs in names(corrs)) {
        args <- c(args, lapply(names(corrs[[cs]]), function(gp) c(cs, gp)))
    }
    # print(args)

    #' Calculate the cost (A * x + e - b)^2 (element-wise) for each case and group
    #'
    #' @param case A string, the case name
    #' @param group A string, the group name (typically 1, 2, 3, ...)
    #' @return A numeric, the cost for the case and group
    ew_cost <- function(case, group) {
        samples <- names(groupmat)[as.character(groupmat[case, ]) == group]
        corrs <- corrs[[case]][[group]]
        out <- 0
        for (reg in rownames(corrs)) {
            if (ncol(corrs) == 1) {
                cors <- corrs[reg, ]
                names(cors) <- colnames(corrs)
            } else {
                cors <- corrs[reg, ]
            }
            cors <- cors[!is.na(cors)]
            targets <- names(cors)
            rexp <- x[samples, reg]
            rexp[rexp < 0] <- 0
            # spread the targets into columns
            texp <- rexp %*% matrix(cors, nrow = 1)
            texp[texp < 0 & !is.na(texp)] <- 1 + texp[texp < 0 & !is.na(texp)]
            cost <- texp - x[samples, targets]
            out <- out + sum(cost^2)
        }
        out
    }

    #' Calculate the cost (cor(A, B) - c0)^2 for each case and group
    #'
    #' @param case A string, the case name
    #' @param group A string, the group name (typically 1, 2, 3, ...)
    #' @return A numeric, the cost for the case and group
    cor_cost <- function(case, group) {
        samples <- names(groupmat)[as.character(groupmat[case, ]) == group]
        corrs <- corrs[[case]][[group]]
        out <- 0
        for (reg in rownames(corrs)) {
            if (ncol(corrs) == 1) {
                cors <- corrs[reg, ]
                names(cors) <- colnames(corrs)
            } else {
                cors <- corrs[reg, ]
            }
            cors <- cors[!is.na(cors)]
            targets <- names(cors)
            rexp <- x[samples, reg]
            # In case there is only one target
            c1 <- correls(rexp, matrix(x[samples, targets], ncol = length(targets)))[, 1]
            out <- out + (c1 - cors)^2
        }
        out
    }

    costfn <- match.arg(costfn, c("ew", "cor"))
    costfn <- ifelse(costfn == "ew", ew_cost, cor_cost)
    sum(unlist(mclapply(args, function(arg) {
        costfn(arg[[1]], arg[[2]])
    })))
    # sum(unlist(lapply(args, function(arg) {
    #     costfn(arg[[1]], arg[[2]])
    # })))
}

#' Simulate dataset
#'
#' @param groupmat A group matrix. Columns are samples and rows are cases.
#'   The values are usually starting from 1 and represent the group of the samples.
#' @param ngenes An integer, total number of genes to simulate
#' @param npair_range An integer vector of length 1 or 2, range of number of gene pairs
#'   for each case
#' @param nregulator_range An vector of length 1 or 2, fraction of all genes to be
#'   regulators
#' @param ntarget_range An integer vector of length 1 or 2, range of number of targets
#'   for each regulator
#' @param gene_prefix A string, prefix for gene names
#' @param gene_index_start An integer, starting index for gene names
#' @param noise A numeric, the standard deviation of the noise to add to the expression
#'   matrix
#' @param maxiter An integer, the maximum number of iterations for optimization
#' @param costfn A string, the cost function to use. It can be "ew" for element-wise
#'   cost, or "cor" for correlation cost.
#' @param ncores An integer, the number of cores to use for parallel computing
#' @param verbose A logical, whether to print the optimization progress
#' @return A list with the following elements:
#'   - emat: A numeric matrix, the expression matrix
#'   - regulator_targets: A data frame with 3 columns. The first column is the
#'     regulator, the second column is empty, and the third column is a string of
#'     targets separated by tabs.
#'   - case_regulators: A data frame with 3 columns. The first column is the case name,
#'     the second column is empty, and the third column is a string of regulators
#'     separated by tabs.
#' @export
sim_dataset <- function(
    groupmat,
    ngenes,
    npair_range = c(10, 16),
    nregulator_range = c(.1, .2),
    ntarget_range = c(10, 25),
    gene_prefix = "Gene",
    gene_index_start = 1,
    noise = 0.05,
    maxiter = 100,
    costfn = "ew",
    ncores = 1,
    verbose = FALSE
) {
    nsamples <- ncol(groupmat)

    genes <- paste0(gene_prefix, (1:ngenes) + gene_index_start - 1)
    ngenes <- length(genes)
    if (length(nregulator_range) == 1) { nregulator_range <- rep(nregulator_range, 2) }
    # Select genes as regulators
    reg_genes <- sample(
        genes,
        size = as.integer(ngenes * runif(1, nregulator_range[1], nregulator_range[2])),
        replace = FALSE
    )
    rest_genes <- genes[!genes %in% reg_genes]
    if (length(ntarget_range) == 1) { ntarget_range <- rep(ntarget_range, 2) }
    # an adjacency matrix for regulator -> targets
    gene_pairs <- matrix(0, nrow = length(reg_genes), ncol = length(rest_genes))
    gene_pairs <- as.data.frame(gene_pairs)
    colnames(gene_pairs) <- rest_genes
    rownames(gene_pairs) <- reg_genes
    for (reg in reg_genes) {
        targets <- sample(
            rest_genes,
            size = sample(ntarget_range[1]:ntarget_range[2], 1),
            replace = FALSE
        )
        gene_pairs[reg, targets] <- 1
    }
    # remove columns with all zeros
    gene_pairs <- gene_pairs[, colSums(gene_pairs) > 0, drop = FALSE]
    # remove rows with all zeros
    gene_pairs <- gene_pairs[rowSums(gene_pairs) > 0, , drop = FALSE]

    if (length(npair_range) == 1) { npair_range <- rep(npair_range, 2) }
    params <- generate_params(groupmat, npair_range, gene_pairs)

    emat <- rbeta(ngenes * nsamples, 2, 2)
    emat[emat < 0] <- 0
    emat[emat > 1] <- 1
    emat <- matrix(emat, ncol = ngenes)
    colnames(emat) <- genes
    rownames(emat) <- colnames(groupmat)
    optgenes <- c(rownames(gene_pairs), colnames(gene_pairs))
    exprs <- emat[, optgenes, drop = FALSE]

    weights <- exprs / exprs  # initial weights
    # x <- cost_fn(weights, exprs, params, groupmat, costfn = costfn, ncores = ncores)
    # print(x)
    # stop()
    opt <- optimg(
        weights, fn = cost_fn,
        x = exprs, corrs = params, groupmat = groupmat,
        costfn = costfn, ncores = ncores,
        method = "ADAM", maxit= maxiter, verbose = verbose
    )
    emat[, optgenes] <- opt$par * exprs

    noise = rnorm(nrow(emat) * ncol(emat), 0, noise)
    noise = matrix(noise, nrow = nrow(emat), byrow = T)
    list(
        emat = emat + noise,
        regulator_targets = genepairs_to_gmt(gene_pairs),
        case_targets = params_to_gmt(params),
        params = params
    )
}

log_info("Reading genotype matrix ...")
groupmat <- read.table(genofile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if (min(groupmat) == 0) { groupmat <- groupmat + 1 }
log_info("Running simulation ...")
x <- sim_dataset(
    groupmat,
    ngenes = ngenes,
    npair_range = npair_range,
    nregulator_range = nregulator_range,
    ntarget_range = ntarget_range,
    gene_prefix = gene_prefix,
    gene_index_start = gene_index_start,
    noise = noise,
    verbose = TRUE
)

emat <- x$emat
if (transpose_exprs) { emat <- t(emat) }

log_info("Saving results ...")
write.table(emat, exfile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(x$case_targets, snpgenefile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(x$regulator_targets, tftargetfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# save the truth
# The truth is in params
# $SNP_19_D
# $SNP_19_D$`1`
#             Gene3     Gene5     Gene15
# Gene42         NA 0.5224003         NA
# Gene44  0.9495537        NA         NA
# Gene43 -0.3960198        NA -0.2115226

# $SNP_19_D$`2`
#            Gene3      Gene5     Gene15
# Gene42        NA -0.7555059         NA
# Gene44 0.6974368         NA         NA
# Gene43 0.3128844         NA -0.1770625

# $SNP_19_D$`3`
#             Gene3      Gene5    Gene15
# Gene42         NA -0.5325191        NA
# Gene44 -0.8119182         NA        NA
# Gene43  0.1202350         NA 0.4867512

# The above params should be converted to
# SNP TF Target Corr1 Corr2 Corr3
# SNP_19_D Gene42 Gene5 0.5224003 -0.7555059 -0.5325191
# SNP_19_D Gene44 Gene3 0.9495537 0.6974368 -0.8119182
# SNP_19_D Gene43 Gene3 -0.3960198 0.3128844 0.1202350
# SNP_19_D Gene43 Gene15 -0.2115226 -0.1770625 0.4867512
truth_data <- function(snp, TF, Target) {
    ugroups <- as.character(1:3)
    do_one_group <- function(g) {
        g <- as.character(g)
        ex_corr_df <- x$params[[snp]][[g]]
        if (is.null(ex_corr_df)) {
            ex_corr <- NA
            ob_corr <- NA
            ob_p <- NA
            plotdf <- data.frame(
                TF = NA,
                Target = NA,
                Group = g,
                cor = ob_corr
            )
        } else {
            ex_corr <- ex_corr_df[TF, Target]
            samples <- names(groupmat)[as.character(groupmat[snp, ]) == g]
            cor <- tryCatch({
                cor.test(x$emat[samples, TF], x$emat[samples, Target])
            }, error = function(e) {
                list(estimate = NA, p.value = NA)
            })
            # cor <- cor.test(x$emat[samples, TF], x$emat[samples, Target])
            ob_corr <- cor$estimate
            ob_p <- cor$p.value
            plotdf <- data.frame(
                TF = x$emat[samples, TF],
                Target = x$emat[samples, Target],
                Group = g,
                cor = ob_corr
            )
        }
        list(
            ex_corr = ex_corr,
            ob_corr = ob_corr,
            ob_p = ob_p,
            plotdf = plotdf)
    }
    g1 <- do_one_group(1)
    g2 <- do_one_group(2)
    g3 <- do_one_group(3)
    plotdf <- rbind(g1$plotdf, g2$plotdf, g3$plotdf)
    max_val <- 1
    # scatter plots faceted by group
    # also add r^2 in the plot
    plotfile <- file.path(truthplotdir, paste0(snp, "-", TF, "-", Target, ".png"))
    p <- ggplot(plotdf, aes(x = TF, y = Target)) +
        geom_point() +
        facet_wrap(~ Group) +
        geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
        geom_text(aes(label = paste0("cor = ", round(cor, 2))), x = 0, y = 0, hjust = 0, vjust = 1) +
        theme_bw() +
        xlim(0, max_val) +
        ylim(0, max_val) +
        labs(title = paste0(snp, "/", TF, "/", Target))
    png(plotfile, res = 70, width = 1400, height = 500)
    print(p)
    dev.off()

    data.frame(
        SNP = snp,
        TF = TF,
        Target = Target,
        ExCorr1 = g1$ex_corr,
        ExCorr2 = g2$ex_corr,
        ExCorr3 = g3$ex_corr,
        ObCorr1 = g1$ob_corr,
        ObCorr2 = g2$ob_corr,
        ObCorr3 = g3$ob_corr,
        ObPval1 = g1$ob_p,
        ObPval2 = g2$ob_p,
        ObPval3 = g3$ob_p,
        AnySig = any(na.omit(c(g1$ob_p, g2$ob_p, g3$ob_p)) < 0.05),
        MaxObCorrDiff = max(abs(
            c(g1$ob_corr, g1$ob_corr, g2$ob_corr) -
            c(g2$ob_corr, g3$ob_corr, g3$ob_corr)
        ))
    )
}

truth <- c()
for (snp in names(x$params)) {
    idx <- which(!is.na(x$params[[snp]][[1]]), arr.ind = TRUE)
    for (i in seq_len(nrow(idx))) {
        TF_idx <- idx[i, 1]
        Target_idx <- idx[i, 2]
        TF = rownames(x$params[[snp]][[1]])[TF_idx]
        Target = colnames(x$params[[snp]][[1]])[Target_idx]
        tmp <- truth_data(snp, TF, Target)
        truth <- rbind(truth, truth_data(snp, TF, Target))
    }
}
truth <- as.data.frame(truth)

write.table(truth, truthallfile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(
    truth[truth$AnySig & truth$MaxObCorrDiff >= 0.5, , drop = FALSE],
    truthfile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
)
