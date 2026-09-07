## Simulation functions for comparing multiple-comparisons corrections
## (uncorrected, Bonferroni-type/Sidak, Holm, single-step max-|T|) across
## model-simplification scenarios (no selection, unrestricted step(),
## interactions-only step()). Reuses step_default()/step_interactions_only()
## from simulate_fig1.R, and maxT_crit() from maxT_correction.R.

## Per-test threshold solving alpha' = 1 - (1-alpha_per_test)^k, the exact
## family-wise rate under independence (the Sidak correction; not literal
## Bonferroni, which instead scales alpha down to alpha/k).
sidak_alpha <- function(alpha, k) 1 - (1 - alpha)^(1 / k)

## Screens `fit`'s terms under all four significance criteria, calibrated
## against `ref` -- the pre-selection full model whose k, correlation
## structure, and residual df define the actual multiple-testing exposure.
## This matters whenever `fit` is itself the result of a selection process
## performed on `ref` (the "selected" scenarios): correcting only for the
## terms that happened to survive selection would ignore the search that
## produced them. Forstmeier & Schielzeth's Fig. 3b makes the same choice,
## applying the full model's Bonferroni threshold to the minimal model
## rather than recomputing it from the minimal model's own term count. When
## `fit` and `ref` are the same object (the "unselected" scenario), this
## reduces to the ordinary single-model correction.
## Returns both the experimentwise indicator (sig_*: is at least one of
## fit's terms significant) and the raw count (n_sig_*: how many of fit's
## terms are significant) for each criterion. A per-parameter error rate
## divides n_sig_* by the *original* full-model k (the `k` column already
## in the condition grid), not by fit's own -- possibly reduced -- term
## count.
screen_criteria <- function(fit, ref = fit, alpha = 0.05, maxT_nsim = 5000) {
  tt <- broom::tidy(fit)
  tt <- tt[tt$term != "(Intercept)", , drop = FALSE]

  ref_tt <- broom::tidy(ref)
  k_ref <- sum(ref_tt$term != "(Intercept)")

  V_ref <- stats::vcov(ref)[-1, -1, drop = FALSE]
  R_ref <- stats::cov2cor(V_ref)
  crit_maxT <- maxT_crit(R_ref, df = ref$df.residual, nsim = maxT_nsim, level = 1 - alpha)

  sig_by_term <- tibble::tibble(
    raw  = tt$p.value < alpha,
    bonf = tt$p.value < sidak_alpha(alpha, k_ref),
    holm = stats::p.adjust(tt$p.value, method = "holm", n = k_ref) < alpha,
    maxT = abs(tt$statistic) > crit_maxT
  )

  ## summarise(across()) rather than pivot_longer(): a fully-pruned model
  ## (tt has 0 rows) would vanish entirely under pivot_longer before
  ## group_by/summarise could aggregate it, silently dropping all four
  ## criterion columns instead of yielding FALSE/0 for each.
  sig_by_term |>
    dplyr::summarise(dplyr::across(
      dplyr::everything(),
      list(sig = ~ any(.x, na.rm = TRUE), n_sig = ~ sum(.x, na.rm = TRUE)),
      .names = "{.fn}_{.col}"
    ))
}

## Metrics unaffected by k_ref (raw p-values aren't corrected), so naive and
## correct calibration give identical values -- kept unsuffixed. The rest
## depend on k_ref and get both a "_naive" and a "_correct" column.
unsuffixed_metrics <- c("sig_raw", "n_sig_raw")
calibrated_metrics <- c("sig_bonf", "sig_holm", "sig_maxT", "n_sig_bonf", "n_sig_holm", "n_sig_maxT")

## Renames and combines one selected scenario's naive/correct screen_criteria()
## outputs into the suffixed column block one_rep_corrections() appends for
## that scenario.
selected_scenario_block <- function(nm, naive, correct) {
  rename_suffix <- function(x, suffix) {
    names(x) <- paste0(names(x), suffix)
    x
  }
  dplyr::bind_cols(
    rename_suffix(correct[unsuffixed_metrics], paste0("_", nm)),
    rename_suffix(naive[calibrated_metrics],   paste0("_", nm, "_naive")),
    rename_suffix(correct[calibrated_metrics], paste0("_", nm, "_correct"))
  )
}

## One simulation replicate: fit the full model once (the "unselected"
## scenario), screen it, then simplify it with each of `simplify_fns` and
## screen those too (calibrated against the full model) -- one column block
## per scenario, suffixed by name. For each selected scenario, corrections
## are computed two ways: "naive" (calibrated on the selected model's own
## surviving term count -- what a user who doesn't think about the search
## that produced it would do) and "correct" (calibrated on the full
## pre-selection model, as in Forstmeier & Schielzeth's Fig. 3b).
##
## At low N relative to k (e.g. N=30-50 with interactions), the
## unorthogonalized random design occasionally produces a near-singular
## design matrix -- not caught by lm()'s aliasing check, but bad enough that
## vcov() picks up non-finite entries and crashes the eigen decomposition
## inside maxT_crit(). Since every draw is null-hypothesis-true regardless,
## a degenerate draw is simply discarded and redrawn.
one_rep_corrections <- function(N, m, interactions, alpha = 0.05, maxT_nsim = 5000,
                                 simplify_fns = list(), max_tries = 20) {
  ## Set here rather than once in the calling script: future::multisession
  ## workers are separate R processes that don't inherit options() set in
  ## the main process (unlike RhpcBLASctl's BLAS thread count under forking,
  ## which is also worker-process-local for the same reason -- see
  ## ~/.claude/r-parallelization.md). Sum-to-zero contrasts make each
  ## main-effect test, in the presence of interactions, a test of the
  ## effect at the population mean of the other factors rather than at an
  ## arbitrary reference level.
  options(contrasts = c("contr.sum", "contr.poly"))
  for (attempt in seq_len(max_tries)) {
    result <- tryCatch({
      dat <- data.frame(y = rnorm(N), make_predictors(N, m))
      ## do.call embeds the data by value in fit_full$call, rather than as a
      ## `dat` symbol: step()'s custom-scope path re-evaluates this call via
      ## update()'s eval(call, parent.frame()), and parent.frame() there is
      ## step()'s own frame, which can't see a `dat` local to this function.
      fit_full <- do.call("lm", list(formula = build_formula(m, interactions), data = dat))

      out <- screen_criteria(fit_full, alpha = alpha, maxT_nsim = maxT_nsim)
      names(out) <- paste0(names(out), "_unselected")

      for (nm in names(simplify_fns)) {
        fit_i <- simplify_fns[[nm]](fit_full, m, interactions)
        naive   <- screen_criteria(fit_i, alpha = alpha, maxT_nsim = maxT_nsim)
        correct <- screen_criteria(fit_i, ref = fit_full, alpha = alpha, maxT_nsim = maxT_nsim)
        out <- dplyr::bind_cols(out, selected_scenario_block(nm, naive, correct))
      }
      out
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
  }
  stop("one_rep_corrections: ", max_tries, " consecutive degenerate design draws for N=",
       N, ", m=", m, ", interactions=", interactions, " -- something is wrong beyond bad luck.")
}

## Run nrep replicates for one (N, m, interactions) condition, in parallel.
run_condition_corrections <- function(N, m, interactions, nrep, simplify_fns = list(),
                                       alpha = 0.05, maxT_nsim = 5000) {
  furrr::future_map_dfr(
    seq_len(nrep),
    function(i) one_rep_corrections(N, m, interactions, alpha, maxT_nsim, simplify_fns),
    .options = furrr::furrr_options(seed = TRUE)
  )
}
