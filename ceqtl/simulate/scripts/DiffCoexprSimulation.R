source("{{biopipen_dir}}/utils/misc.R")

library(foreach)
library(doParallel)
library(nleqslv)

genofile <- {{in.genofile | r}}
exfile <- {{out.exfile | r}}
snpgenefile <- {{out.snpgenefile | r}}
tftargetfile <- {{out.tftargetfile | r}}

seed <- {{envs.seed | r}}
ncores <- {{envs.ncores | r}}
ngenes <- {{envs.ngenes | r}}
npair_range <- {{envs.npair_range | r}}
nregulator_range <- {{envs.nregulator_range | r}}
ntarget_range <- {{envs.ntarget_range | r}}
gene_prefix <- {{envs.gene_prefix | r}}
gene_index_start <- {{envs.gene_index_start | r}}
bionoise <- {{envs.bionoise | r}}
expnoise <- {{envs.expnoise | r}}
transpose_exprs <- {{envs.transpose_exprs | r}}

set.seed(seed)
registerDoParallel(cores=24)

#### Generate genotype matrix
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
        probs <- 0
        while (min(probs) < minprob) { probs <- rnorm(ngroups) }
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
            while (min(probs) < minprob) { probs <- probs_fn(ngroups) }
        } else if (length(probs) == ngroups - 1) {
            probs <- c(probs, 1 - sum(probs))
        }
        groups <- sample(1:ngroups, size = nsamples, prob = probs, replace = TRUE)
        groupmat <- rbind(groupmat, groups)
        rownames(groupmat)[nrow(groupmat)] = n
    }
    colnames(groupmat) <- samples

    return(groupmat)
}

as_groupmat.data.frame <- function(x, ...) { x }

invisible(setGeneric("as_groupmat", function(x, ...) standardGeneric("as_groupmat")))
setMethod("as_groupmat", "list", as_groupmat.list)
setMethod("as_groupmat", "integer", as_groupmat.int)
setMethod("as_groupmat", "numeric", as_groupmat.int)
setMethod("as_groupmat", "data.frame", as_groupmat.data.frame)
### End of as_groupmat



add_bionoise <- function(x, lnorm){
    if (all(lnorm == 1)) {
        return(x)
    }

    #transformation fn
    f <- function(a, b){
        fnval = a * exp(-0.01 *( 1 - a/(1 - a))) - b
        return(fnval)
    }
    optf <- function(a, b){
        fnval = a * exp(-0.01 *( 1 - a/(1 - a))) - b
        return(abs(fnval))
    }

    #transformation
    y = x
    tfx = x>=0.65
    y[tfx] = f(x[tfx], 0)
    y[!tfx] = x[!tfx]

    #apply lognormal noise
    newy = y * lnorm

    #inverse transformation
    newx = x
    tfy = newy>=0.65
    newx[tfy] = unlist(sapply(newy[tfy], function (a)
        optim(
        runif(1),
        optf,
        b = a,
        lower = 0,
        upper = 1,
        method = 'Brent',
        control = list('abstol' = 1E-8)
        )$par))
    newx[!tfy] = newy[!tfy]
    return(newx)
}

corr_fun <- function(x, slope, intercept){
    return(slope * x + intercept)
}

ode <- function(exprs, slopes, intercepts, lnnoise){
    rates <- exprs
    pairs <- names(slopes)
    for (pair in pairs) {
        gs <- strsplit(pair, " -> ")[[1]]
        slope <- slopes[pair]
        intercept <- intercepts[pair]
        if (slope > 0) {
            rates[gs[2]] <- add_bionoise(
                corr_fun(exprs[gs[1]], slope, intercept), lnnoise[gs[2]]
            ) - rates[gs[2]]
        } else {
            rates[gs[2]] <- add_bionoise(
                1 - corr_fun(exprs[gs[1]], -slope, intercept), lnnoise[gs[2]]
            ) - rates[gs[2]]
        }
    }
    return(rates)
}

generate_params <- function(groupmat, npair_range, gene_pairs) {
    out = list()
    # list(
    #    c1 = list( "g1 -> g2": list(1 => c(slope, intercept), 2 => c(slope, intercept)) )
    # )
    for (cs in rownames(groupmat)) {
        ugroups <- unique(unlist(groupmat[cs, ]))
        ngroups <- length(ugroups)
        npairs <- sample(npair_range[1]:npair_range[2], 1)
        gps <- sample(gene_pairs, npairs, replace = FALSE)
        out[[cs]] = list()
        for (gp in gps) {
            slopes <- runif(ngroups, -1, 1)
            while(max(apply(combn(slopes, 2), 2, function(x) abs(x[1] - x[2]))) < 0.5) {
                slopes <- runif(ngroups, -1, 1)
            }
            intersects <- runif(ngroups, -1, 1)
            out[[cs]][[gp]] = lapply(seq_len(ngroups), function(i) c(slopes[i], intersects[i]))
        }
    }

    return(out)
}

get_params_by_sample <- function(pms, samplegroup) {
    ## params
    # list( "g1 -> g2": list(1 => c(slope, intercept), 2 => c(slope, intercept)) )
    pairs = names(pms)
    slopes = sapply(pairs, function(x) pms[[x]][[samplegroup]][1])
    intercepts = sapply(pairs, function(x) pms[[x]][[samplegroup]][2])
    names(slopes) = pairs
    names(intercepts) = pairs
    return(list(slopes = slopes, intercepts = intercepts))
}

genepairs_to_gmt <- function(gene_pairs) {
    # convert to gmt dataframe
    # gene_pairs = c("g1 -> g2", "g1 -> g3", "g2 -> g3")
    out <- list()
    for (pair in gene_pairs) {
        genes <- strsplit(pair, " -> ")[[1]]
        out[[genes[1]]] <- c(out[[genes[1]]], genes[2])
    }
    return (data.frame(
        V1 = names(out),
        V2 = "",
        V3 = unlist(lapply(out, function(x) paste(x, collapse = "\t")))
    ))
}

params_to_gmt <- function(params) {
    # convert to gmt dataframe
    # params = list(
    #    c1 = list( "g1 -> g2": list(1 => c(slope, intercept), 2 => c(slope, intercept)) )
    # )
    out <- list()
    for (cs in names(params)) {
        pairs <- names(params[[cs]])
        regulators <- unique(sapply(pairs, function(x) strsplit(x, " -> ")[[1]][1]))
        out[[cs]] <- paste(regulators, collapse = "\t")
    }
    return (data.frame(
        V1 = names(out),
        V2 = "",
        V3 = unlist(out)
    ))
}

sim_dataset <- function(
    groupmat,
    ngenes,
    npair_range = c(10, 16),
    nregulator_range = c(.1, .2),
    ntarget_range = c(10, 25),
    gene_prefix = "Gene",
    gene_index_start = 1,
    bionoise = 0.05,
    expnoise = 0.05
) {
    nsamples <- ncol(groupmat)

    genes <- paste0(gene_prefix, (1:ngenes) + gene_index_start - 1)
    ngenes <- length(genes)
    if (length(nregulator_range) == 1) { nregulator_range <- rep(nregulator_range, 2) }
    reg_genes <- sample(
        genes,
        size = as.integer(ngenes * runif(1, nregulator_range[1], nregulator_range[2])),
        replace = FALSE
    )
    if (length(ntarget_range) == 1) { ntarget_range <- rep(ntarget_range, 2) }

    log_info("Generating regulator-gene pairs ...")
    gene_pairs <- unlist(lapply(reg_genes, function(reg_gene) {
        targets <- sample(
            genes,
            size = sample(ntarget_range[1]:ntarget_range[2], 1),
            replace = FALSE
        )
        paste0(reg_gene, " -> ", targets)
    }))

    if (length(npair_range) == 1) { npair_range <- rep(npair_range, 2) }

    log_info("Randomizing parameters for gene pairs ...")
    params <- generate_params(groupmat, npair_range, gene_pairs)

    lnnoise = exp(rnorm(nsamples * ngenes, 0, bionoise))
    lnnoise = matrix(lnnoise, nrow = nsamples, byrow = T)
    colnames(lnnoise) = genes

    log_info("Initializing gene expression ...")
    emat <- rbeta(ngenes * nsamples, 2, 2)
    emat[emat < 0] <- 0
    emat[emat > 1] <- 1
    emat <- matrix(emat, nrow = ngenes)
    rownames(emat) <- genes

    log_info("Adjusting gene expression for differential coexpression ...")
    for (cs in names(params)) {
        log_info("- case: {cs}")
        emat = foreach(i = seq_len(nsamples), .packages = c('nleqslv'), .combine = cbind) %dopar% {
        # tmp <- lapply(seq_len(nsamples), function(i) {
            # print(i)
            parms <- get_params_by_sample(params[[cs]], groupmat[cs, i])
            soln = suppressWarnings(nleqslv(
                emat[, i],
                ode,
                slopes = parms$slopes, intercepts = parms$intercepts, lnnoise = lnnoise[i, ]
            ))
            # return(c(soln$x, soln$termcd))
            return(soln$x)
        }
        # print(length(tmp))
        # emat = do.call(cbind, tmp)
        # print(dim(emat))
        # # emat = cbind(emat)
        # print(emat)
    }

    # termcd = res[nrow(res),]
    # emat = res[-(nrow(res)), , drop = F]
    colnames(emat) = colnames(groupmat)
    rownames(emat) = genes

    # if (!all(termcd == 1)) {
    #     nc = termcd != 1
    #     msg = 'Simulations for the following samples did not converge:'
    #     sampleids = paste(colnames(emat)[nc], ' (', termcd[nc], ')', sep = '', collapse = ', ')
    #     msg = paste(c(msg, sampleids), collapse = '\n\t')
    #     msg = paste(msg, 'format: sampleid (termination condition)', sep = '\n\n\t')
    #     warning(msg)

    #     emat = emat[, !nc]
    # }

    expnoise = rnorm(nrow(emat) * ncol(emat), 0, expnoise)
    expnoise = matrix(expnoise, nrow = nrow(emat), byrow = T)
    list(
        emat = emat + expnoise,
        regulator_targets = genepairs_to_gmt(gene_pairs),
        case_regulators = params_to_gmt(params)
    )
}

log_info("Reading genotype matrix ...")
groupmat <- read.table(genofile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if (min(groupmat) == 0) { groupmat <- groupmat + 1 }
x <- sim_dataset(
    groupmat,
    ngenes = ngenes,
    npair_range = npair_range,
    nregulator_range = nregulator_range,
    ntarget_range = ntarget_range,
    gene_prefix = gene_prefix,
    gene_index_start = gene_index_start,
    bionoise = bionoise,
    expnoise = expnoise
)

emat <- x$emat
if (transpose_exprs) { emat <- t(emat) }

log_info("Saving results ...")
write.table(emat, exfile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(x$case_regulators, snpgenefile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(x$regulator_targets, tftargetfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
