## Benchmark scaling of phylogenetic covariance/precision matrix construction
## (vcv/solve vs MRFtools::mrf_penalty + sparsification) with subsample size,
## using the full fishtree phylogeny as the source tree.

library(fishtree)
library(ape)
library(MRFtools)
library(Matrix)
library(future)
library(furrr)
library(here)

## avoid BLAS oversubscription when forking workers (each fork would
## otherwise also try to multithread its own BLAS calls)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
}

set.seed(20260817)

n_workers <- 10
plan(multicore, workers = n_workers)

sizes <- c(100, 500, 1000, 2500, 5000, 7500, 10000)
n_reps <- 10
outfile <- here("phyloslopes", "phyloslopes_testinv.rds")

f_full <- fishtree_phylogeny()
message(sprintf("full tree: %d tips, %d internal nodes", length(f_full$tip.label), f_full$Nnode))

time_it <- function(expr) {
    t0 <- proc.time()
    val <- expr
    dt <- proc.time() - t0
    list(value = val,
         user = unname(dt[["user.self"]]),
         sys = unname(dt[["sys.self"]]),
         elapsed = unname(dt[["elapsed"]]))
}

step_row <- function(step, timed, size_bytes, sparsity = NA_real_) {
    data.frame(step = step,
               user = timed$user, sys = timed$sys, elapsed = timed$elapsed,
               size_bytes = size_bytes, sparsity = sparsity)
}

## run all 6 timed steps for one random subsample of `n_tips` tips
run_one <- function(n_tips, rep_id, full_tree) {
    tips_sub <- sample(full_tree$tip.label, n_tips)
    tr <- keep.tip(full_tree, tips_sub)

    rows <- list()

    r <- time_it(vcv(tr))
    v <- r$value
    rows[["vcv"]] <- step_row("vcv", r, as.numeric(object.size(v)))

    r <- time_it(solve(v))
    rows[["solve"]] <- step_row("solve", r, as.numeric(object.size(r$value)))
    rm(v)

    r <- time_it(mrf_penalty(tr, internal_nodes = FALSE))
    p_tips <- r$value
    rows[["mrf_tips"]] <- step_row("mrf_tips", r, as.numeric(object.size(p_tips)))

    r <- time_it(Matrix(p_tips, sparse = TRUE))
    p_tips_s <- r$value
    rows[["mrf_tips_sparse"]] <- step_row("mrf_tips_sparse", r,
                                           as.numeric(object.size(p_tips_s)),
                                           sparsity = mean(p_tips_s != 0))
    rm(p_tips, p_tips_s)

    r <- time_it(mrf_penalty(tr, internal_nodes = TRUE))
    p_full <- r$value
    rows[["mrf_full"]] <- step_row("mrf_full", r, as.numeric(object.size(p_full)))

    r <- time_it(Matrix(p_full, sparse = TRUE))
    p_full_s <- r$value
    rows[["mrf_full_sparse"]] <- step_row("mrf_full_sparse", r,
                                           as.numeric(object.size(p_full_s)),
                                           sparsity = mean(p_full_s != 0))
    rm(p_full, p_full_s)

    out <- do.call(rbind, rows)
    out$n_tips <- n_tips
    out$rep <- rep_id
    out$actual_ntip <- tr$Nnode + 1L
    out$actual_nnode <- tr$Nnode
    rownames(out) <- NULL
    out[c("n_tips", "rep", "actual_ntip", "actual_nnode", "step",
          "user", "sys", "elapsed", "size_bytes", "sparsity")]
}

run_one_safe <- function(n_tips, rep_id, full_tree) {
    tryCatch(
        run_one(n_tips, rep_id, full_tree),
        error = function(e) {
            data.frame(n_tips = n_tips, rep = rep_id,
                       actual_ntip = NA_integer_, actual_nnode = NA_integer_,
                       step = "ERROR", user = NA_real_, sys = NA_real_,
                       elapsed = NA_real_, size_bytes = NA_real_,
                       sparsity = NA_real_,
                       error_message = conditionMessage(e))
        }
    )
}

## resume support: skip sizes already checkpointed in outfile
results <- if (file.exists(outfile)) readRDS(outfile) else list()
done_sizes <- suppressWarnings(as.integer(names(results)))

for (n_tips in sizes) {
    if (n_tips %in% done_sizes) {
        message(sprintf("size %d already done, skipping", n_tips))
        next
    }
    message(sprintf("size %d: running %d reps across %d workers...", n_tips, n_reps, n_workers))
    t_size <- system.time({
        reps <- future_map(
            seq_len(n_reps),
            run_one_safe,
            n_tips = n_tips, full_tree = f_full,
            .options = furrr_options(seed = TRUE)
        )
    })
    results[[as.character(n_tips)]] <- do.call(rbind, reps)
    saveRDS(results, outfile)
    message(sprintf("size %d done in %.1f s (checkpointed to %s)", n_tips, t_size[["elapsed"]], outfile))
}

message("all sizes complete")
