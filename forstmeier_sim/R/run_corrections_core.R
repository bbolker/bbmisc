## Shared "run the multiple-comparisons-correction simulation grid" logic
## -- three model-simplification scenarios (unselected full model,
## unrestricted step(), interactions-only step()), each screened under all
## four significance criteria (uncorrected, Dunn-Sidak, Holm, single-step
## max-|T|; the two selected scenarios additionally screened under both
## naive and correct k-counting for the corrections). Parametrized by
## output path so a rerun can be pointed at a scratch file without
## clobbering existing results. Reused by run_corrections_simulation.R.

run_corrections_simulation <- function(results_path, nrep = 5000, n_workers = 24) {
  if (file.exists(results_path)) {
    message(basename(results_path), " already exists -- skipping simulation rerun.")
    return(invisible(NULL))
  }

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

  results <- conditions |>
    dplyr::mutate(
      k = purrr::map2_dbl(m, interactions, n_predictors),
      sims = purrr::pmap(list(N, m, interactions), function(N, m, interactions) {
        run_condition_corrections(N, m, interactions, nrep, simplify_fns = simplify_fns)
      })
    ) |>
    tidyr::unnest(sims)

  future::plan(future::sequential)

  saveRDS(results, results_path)
}
