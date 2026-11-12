## Shared "run the multiple-comparisons-correction simulation grid" logic
## -- three model-simplification scenarios (unselected full model,
## unrestricted step(), interactions-only step()), each screened under all
## five significance criteria (uncorrected, Bonferroni, Dunn-Sidak, Holm,
## single-step max-|T|; the two selected scenarios additionally screened under both
## naive and correct k-counting for the corrections). Parametrized by
## output path so a rerun can be pointed at a scratch file without
## clobbering existing results. Reused by run_corrections_simulation.R.
##
## Each run appends its timing (start/end, elapsed minutes, host, worker
## and core counts, nrep) to corrections_run_log.csv next to results_path,
## and stores the same record as the "run_info" attribute of the saved
## results. Progress is written to corrections_progress.csv (also next to
## results_path, overwritten at the start of each run): one line per
## condition as it finishes, with its elapsed time.
##
## Always recomputes and overwrites results_path. The targets pipeline
## decides when a rerun is needed (whenever any of the simulation code
## files changes); skipping here when the file already exists would let
## stale results survive code changes.

## at most 24 workers, leaving at least one core free (detectCores() can return NA)
run_corrections_simulation <- function(results_path, nrep = 5000,
                                       n_workers = min(24, max(1, parallel::detectCores() - 1, na.rm = TRUE))) {
  start_time <- Sys.time()
  set.seed(20260906)

  ## Single-threaded BLAS/OpenMP inside each worker -- see
  ## ~/.claude/r-parallelization.md.
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)

  ## Cost is dominated by the per-replicate maxT_crit() Monte Carlo draws
  ## (one per scenario for "unselected", two -- naive and correct -- for
  ## each selected scenario) at large k; 24 workers keeps the full grid (3
  ## N values x 6 m x 2 interaction settings) practical.
  future::plan(future::multisession, workers = n_workers)

  simplify_fns <- list(step = step_default, step_restricted = step_interactions_only)

  conditions <- tidyr::expand_grid(
    N = sim_N_values,
    m = 1:6,
    interactions = c(FALSE, TRUE)
  )

  progress_path <- file.path(dirname(results_path), "corrections_progress.csv")
  writeLines("condition,n_conditions,N,m,interactions,finished,elapsed_min", progress_path)

  results <- conditions |>
    dplyr::mutate(
      k = purrr::map2_dbl(m, interactions, n_predictors),
      sims = purrr::pmap(list(seq_along(N), N, m, interactions), function(i, N, m, interactions) {
        t0 <- Sys.time()
        res <- run_condition_corrections(N, m, interactions, nrep, simplify_fns = simplify_fns)
        cat(sprintf("%d,%d,%d,%d,%s,%s,%.2f\n", i, nrow(conditions), N, m, interactions,
                    format(Sys.time(), "%H:%M:%S"),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))),
            file = progress_path, append = TRUE)
        res
      })
    ) |>
    tidyr::unnest(sims)

  future::plan(future::sequential)

  end_time <- Sys.time()
  run_info <- data.frame(
    start = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
    end = format(end_time, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_min = round(as.numeric(difftime(end_time, start_time, units = "mins")), 1),
    host = Sys.info()[["nodename"]],
    n_workers = n_workers,
    n_cores = parallel::detectCores(),
    nrep = nrep
  )
  attr(results, "run_info") <- run_info
  saveRDS(results, results_path)

  log_path <- file.path(dirname(results_path), "corrections_run_log.csv")
  utils::write.table(run_info, log_path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(log_path), append = file.exists(log_path))
  message(sprintf("simulation took %.1f min on %s with %d workers",
                  run_info$elapsed_min, run_info$host, n_workers))
}
