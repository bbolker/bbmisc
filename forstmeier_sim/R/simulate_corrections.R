## Simulation functions for comparing multiple-comparisons corrections
## (uncorrected, Bonferroni-type/Sidak, Holm, single-step max-|T|) across
## model-simplification scenarios (no selection, unrestricted step(),
## interactions-only step()). Reuses step_default()/step_interactions_only()
## from simulate_fig1.R, and maxT_crit() from maxT_correction.R.

## Per-test threshold solving alpha' = 1 - (1-alpha_per_test)^k, the exact
## family-wise rate under independence (the Sidak correction; not literal
## Bonferroni, which instead scales alpha down to alpha/k).
sidak_alpha <- function(alpha, k) 1 - (1 - alpha)^(1 / k)

## TRUE for a main-effect term (no ":" in its name), FALSE for an
## interaction term.
is_main_effect_term <- function(term) !grepl(":", term, fixed = TRUE)

## Tidies `fit`, drops the intercept row, and optionally further restricts
## to the terms `term_subset` selects (a predicate on term names, e.g.
## is_main_effect_term) -- shared by compute_calibration(), screen_criteria()
## and screen_raw().
tidy_fit <- function(fit, term_subset = NULL) {
  tt <- broom::tidy(fit) |> subset(term != "(Intercept)")
  if (!is.null(term_subset)) tt <- tt[term_subset(tt$term), , drop = FALSE]
  tt
}

## Computes the multiple-testing exposure (k, and the single-step max-|T|
## critical value) implied by `ref` -- the pre-selection full model whose
## k, correlation structure, and residual df define the actual exposure,
## regardless of which/how many of its terms are later asked about. This
## matters whenever the model being screened is itself the result of a
## selection process performed on `ref` (the "selected" scenarios):
## correcting only for the terms that happened to survive selection would
## ignore the search that produced them. Forstmeier & Schielzeth's Fig. 3b
## makes the same choice, applying the full model's Bonferroni threshold to
## the minimal model rather than recomputing it from the minimal model's
## own term count.
##
## `ref_subset`, if given, restricts k/the correlation matrix to just the
## subset of ref's own terms it selects (e.g. is_main_effect_term) -- used
## for the "naive" main-effects-only calibration, where a user who only
## cares about main effects would count only those as their exposure.
## Leaving it NULL (the default) uses ref's full term set, as for the
## ordinary and "correct"/maximal calibrations.
compute_calibration <- function(ref, ref_subset = NULL, alpha = 0.05, maxT_nsim = 5000) {
  ref_terms <- tidy_fit(ref, ref_subset)$term
  k_ref <- length(ref_terms)

  V_ref <- stats::vcov(ref)[ref_terms, ref_terms, drop = FALSE]
  R_ref <- stats::cov2cor(V_ref)
  crit_maxT <- maxT_crit(R_ref, df = ref$df.residual, nsim = maxT_nsim, level = 1 - alpha)

  list(k_ref = k_ref, crit_maxT = crit_maxT)
}

## Screens `fit`'s terms (optionally restricted to the subset `test_subset`
## selects, e.g. is_main_effect_term) under all four significance criteria,
## using a `calibration` (see compute_calibration()) computed separately --
## letting the same calibration be reused across multiple test_subsets
## (e.g. all terms and main-effects-only both calibrated on the full
## pre-selection model) without recomputing its expensive max-|T| Monte
## Carlo critical value each time.
##
## Returns both the experimentwise indicator (sig_*: is at least one of
## fit's screened terms significant) and the raw count (n_sig_*: how many
## are). A per-parameter error rate divides n_sig_* by the relevant
## original full-model term count (`k`, or `m` for the main-effects-only
## variants) rather than by fit's own -- possibly reduced -- term count.
screen_criteria <- function(fit, calibration, alpha = 0.05, test_subset = NULL) {
  tt <- tidy_fit(fit, test_subset)
  k_ref <- calibration$k_ref
  crit_maxT <- calibration$crit_maxT

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

## Screens fit's terms (optionally restricted by `subset`) under the
## uncorrected criterion only -- no calibration object needed, since raw
## p-values don't depend on k. Used for the lightweight main-effects-only
## and cross-check metrics, which don't need max-|T|'s Monte Carlo cost.
screen_raw <- function(fit, alpha = 0.05, subset = NULL) {
  tt <- tidy_fit(fit, subset)
  sig <- tt$p.value < alpha
  tibble::tibble(sig_raw = any(sig, na.rm = TRUE), n_sig_raw = sum(sig, na.rm = TRUE))
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
## pre-selection model, as in Forstmeier & Schielzeth's Fig. 3b). The full
## pre-selection model's calibration (calib_full) is computed once and
## reused for "unselected" and every scenario's "correct" -- these would
## otherwise each independently re-estimate the identical max-|T| critical
## value via its own noisy 5000-draw Monte Carlo call.
##
## When `interactions` is TRUE, also tracks the same breakdown restricted
## to just the m main-effect terms (excluding interactions), suffixed
## "_main_<scenario>". Both calibrations change accordingly: "correct_main"
## uses calib_full_main -- fit_full's own calibration, but with k = m and
## the correlation matrix restricted to just its main-effect sub-block,
## since the main-effects-only metrics were never intended to also test
## the interaction terms, so m (not m + choose(m,2)) is the actual
## multiple-testing exposure being corrected for. "naive_main" gets its
## own calibration too, from the selected model's own *surviving* main
## effects -- always m for step_restricted (main effects are never
## dropped), potentially fewer for step. The unselected scenario
## additionally gets an uncorrected-only cross-check: an independently
## fit main-effects-only model (no interaction terms at all), for
## comparing against the main-effect subset of the full interactive
## model's own screening -- these needn't agree replicate-by-replicate
## (the random factor assignment isn't exactly orthogonal in any single
## finite draw) but should match in aggregate over many replicates.
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
  ## effect at the population mean of the other factors rather than an
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

      calib_full <- compute_calibration(fit_full, alpha = alpha, maxT_nsim = maxT_nsim)
      out <- screen_criteria(fit_full, calib_full, alpha = alpha)
      names(out) <- paste0(names(out), "_unselected")

      if (interactions) {
        ## Calibrated on fit_full restricted to just its m main-effect
        ## terms -- k = m, not m + choose(m,2) -- since the main-effects-only
        ## metrics were never intended to also test the interaction terms;
        ## the correlation structure and residual df still come from
        ## fit_full itself (same fit, just querying a subset of its
        ## coefficients), unlike calib_full's own unrestricted k.
        calib_full_main <- compute_calibration(fit_full, ref_subset = is_main_effect_term,
                                                alpha = alpha, maxT_nsim = maxT_nsim)

        main_unselected <- screen_raw(fit_full, alpha = alpha, subset = is_main_effect_term)
        names(main_unselected) <- paste0(names(main_unselected), "_main_unselected")

        fit_mainonly <- do.call("lm", list(formula = build_formula(m, FALSE), data = dat))
        crosscheck <- screen_raw(fit_mainonly, alpha = alpha)
        names(crosscheck) <- paste0(names(crosscheck), "_mainonly_crosscheck")

        out <- dplyr::bind_cols(out, main_unselected, crosscheck)
      }

      for (nm in names(simplify_fns)) {
        fit_i <- simplify_fns[[nm]](fit_full, m, interactions)
        calib_naive <- compute_calibration(fit_i, alpha = alpha, maxT_nsim = maxT_nsim)
        naive   <- screen_criteria(fit_i, calib_naive, alpha = alpha)
        correct <- screen_criteria(fit_i, calib_full, alpha = alpha)
        out <- dplyr::bind_cols(out, selected_scenario_block(nm, naive, correct))

        if (interactions) {
          calib_naive_main <- compute_calibration(fit_i, ref_subset = is_main_effect_term,
                                                   alpha = alpha, maxT_nsim = maxT_nsim)
          naive_main   <- screen_criteria(fit_i, calib_naive_main, alpha = alpha,
                                           test_subset = is_main_effect_term)
          correct_main <- screen_criteria(fit_i, calib_full_main, alpha = alpha,
                                           test_subset = is_main_effect_term)
          out <- dplyr::bind_cols(out, selected_scenario_block(paste0("main_", nm), naive_main, correct_main))
        }
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
