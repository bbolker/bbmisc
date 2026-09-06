## Simulation functions for comparing multiple-comparisons corrections
## (uncorrected, Bonferroni-type/Sidak, Holm, single-step max-|T|) across
## model-simplification scenarios (no selection, unrestricted step(),
## interactions-only step()). Reuses step_default()/step_interactions_only()
## and any_significant() from simulate_fig1.R, and maxT_crit() from
## maxT_correction.R.

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
screen_criteria <- function(fit, ref = fit, alpha = 0.05, maxT_nsim = 5000) {
  tt <- broom::tidy(fit)
  tt <- tt[tt$term != "(Intercept)", , drop = FALSE]

  ref_tt <- broom::tidy(ref)
  k_ref <- sum(ref_tt$term != "(Intercept)")

  V_ref <- stats::vcov(ref)[-1, -1, drop = FALSE]
  R_ref <- stats::cov2cor(V_ref)
  crit_maxT <- maxT_crit(R_ref, df = ref$df.residual, nsim = maxT_nsim, level = 1 - alpha)

  tibble::tibble(
    sig_raw  = any_significant(fit, alpha),
    sig_bonf = any(tt$p.value < sidak_alpha(alpha, k_ref), na.rm = TRUE),
    sig_holm = any(stats::p.adjust(tt$p.value, method = "holm", n = k_ref) < alpha, na.rm = TRUE),
    sig_maxT = any(abs(tt$statistic) > crit_maxT, na.rm = TRUE)
  )
}

## One simulation replicate: fit the full model once (the "unselected"
## scenario), screen it, then simplify it with each of `simplify_fns` and
## screen those too (calibrated against the full model) -- one column block
## per scenario, suffixed by name.
##
## At low N relative to k (e.g. N=30-50 with interactions), the
## unorthogonalized random design occasionally produces a near-singular
## design matrix -- not caught by lm()'s aliasing check, but bad enough that
## vcov() picks up non-finite entries and crashes the eigen decomposition
## inside maxT_crit(). Since every draw is null-hypothesis-true regardless,
## a degenerate draw is simply discarded and redrawn.
one_rep_corrections <- function(N, m, interactions, alpha = 0.05, maxT_nsim = 5000,
                                 simplify_fns = list(), max_tries = 20) {
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

      ## For each selected scenario, compute corrections two ways: "naive"
      ## (calibrated on the selected model's own surviving term count --
      ## what a user who doesn't think about the search that produced it
      ## would do) and "correct" (calibrated on the full pre-selection
      ## model, as in Forstmeier & Schielzeth's Fig. 3b). sig_raw doesn't
      ## depend on k, so it's identical either way and kept unsuffixed.
      for (nm in names(simplify_fns)) {
        fit_i <- simplify_fns[[nm]](fit_full, m, interactions)
        naive   <- screen_criteria(fit_i, alpha = alpha, maxT_nsim = maxT_nsim)
        correct <- screen_criteria(fit_i, ref = fit_full, alpha = alpha, maxT_nsim = maxT_nsim)

        block <- tibble::tibble(
          !!paste0("sig_raw_", nm)             := correct$sig_raw,
          !!paste0("sig_bonf_", nm, "_naive")   := naive$sig_bonf,
          !!paste0("sig_bonf_", nm, "_correct") := correct$sig_bonf,
          !!paste0("sig_holm_", nm, "_naive")   := naive$sig_holm,
          !!paste0("sig_holm_", nm, "_correct") := correct$sig_holm,
          !!paste0("sig_maxT_", nm, "_naive")   := naive$sig_maxT,
          !!paste0("sig_maxT_", nm, "_correct") := correct$sig_maxT
        )
        out <- dplyr::bind_cols(out, block)
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
