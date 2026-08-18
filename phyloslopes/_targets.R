library(targets)

tar_option_set(
  packages = c("ade4", "ape", "dplyr", "glmmTMB", "RTMB", "reformulas",
               "Matrix", "MRFtools", "phyr", "mgcv", "microbenchmark",
               "MASS", "tibble", "ggplot2", "cowplot", "broom.mixed",
               "performance")
)

## Every script in this project (phyloslopes_tests.R, phyloslopes_tiny.R,
## ...) is a plain top-to-bottom script, not a function -- see README.md --
## so each target below just source()s a script for its side effect (it
## writes its own output file(s)) and tracks those outputs with
## format = "file", rather than trying to pass R objects between targets
## directly. Two reasons for that, not just matching the existing style:
##   (1) TMB/RTMB ADFun objects hold C++ pointers and don't survive being
##       cached/reloaded across a session boundary (see
##       phyloslopes_tiny_predcovs.R's comments on this) -- so a target
##       can never safely *return* a fitted model object anyway.
##   (2) these scripts weren't written as functions, and turning them into
##       one would be a much bigger, riskier change than wiring them up as
##       they are.
##
## "Source-code" targets (utils_r, tiny_r, tiny_predcovs_r, linear_r,
## tests_r, bench1_r, bench1_plot_r, qmd_source) just track a file's hash
## so that editing it correctly invalidates every downstream target that
## sources/renders it -- targets' static dependency analysis can't see
## *inside* a source()d script's own nested source()/render() calls, so
## without these, an edit to e.g. phyloslopes_utils.R or phyloslopes_tiny.R
## wouldn't trigger any reruns at all. Every script that gets source()d
## below has one of these, forced as a dependency in that script's own
## target -- the one exception is phyloslopes_testinv.R/_plot.R, which
## aren't wired into this pipeline at all (see the qmd_html comment).
list(
  tar_target(utils_r, "phyloslopes_utils.R", format = "file"),

  ## -- tiny (5-tip) worked-example pipeline ----------------------------
  tar_target(tiny_r, "phyloslopes_tiny.R", format = "file"),
  tar_target(
    phyloslopes_tiny_rda,
    { utils_r; tiny_r; source("phyloslopes_tiny.R"); "phyloslopes_tiny.rda" },
    format = "file"
  ),
  tar_target(tiny_predcovs_r, "phyloslopes_tiny_predcovs.R", format = "file"),
  tar_target(
    phyloslopes_tiny_predcovs_rda,
    {
      utils_r; tiny_predcovs_r
      phyloslopes_tiny_rda  ## predcovs load()s this .rda
      source("phyloslopes_tiny_predcovs.R")
      "phyloslopes_tiny_predcovs.rda"
    },
    format = "file"
  ),

  ## -- linear (log_bm) models: null/additive/independent/separable ----
  ## null/additive/independent are plain glmmTMB fits; separable reuses
  ## the RTMB nllfun_sep machinery (propto() can't express its joint
  ## phylo x 2x2-trait Kronecker covariance -- see phyloslopes_linear.R's
  ## header). Not expensive (a handful of fits on 49 data points), so no
  ## cue = "never" needed. `independent` doesn't converge cleanly on this
  ## un-replicated dataset (documented in the script) -- kept anyway.
  tar_target(linear_r, "phyloslopes_linear.R", format = "file"),
  tar_target(
    phyloslopes_linear_rda,
    { utils_r; linear_r; source("phyloslopes_linear.R"); "phyloslopes_linear.rda" },
    format = "file"
  ),

  ## -- spline (log_bm) models: null/additive/separable/tensor ----------
  ## analogous to the linear-models block above but for a smooth (not
  ## just linear) log_bm trend -- null is a plain glmmTMB fit; additive/
  ## separable/tensor are all RTMB (nllfun_spline_*), even though additive
  ## alone could be glmmTMB (phyloslopes.qmd's fit_spline_glmmTMB
  ## validates this), so all three share the same b_spline/b_phylo
  ## machinery and stay directly nested/comparable. See
  ## phyloslopes_combo.R's header for the full rationale -- most of its
  ## setup/fitting code is adapted directly from phyloslopes.qmd's
  ## "phylogenetic splines" section.
  tar_target(combo_r, "phyloslopes_combo.R", format = "file"),
  tar_target(
    phyloslopes_combo_rda,
    { utils_r; combo_r; source("phyloslopes_combo.R"); "phyloslopes_combo.rda" },
    format = "file"
  ),

  ## -- 49-species consistency checks + benchmarks ----------------------
  ## EXPENSIVE (~15 model refits + consistency checks, ~10+ min): cue =
  ## "never" so it isn't rerun on every tar_make() -- but targets still
  ## builds a cue = "never" target the *first* time (no recorded output
  ## yet); after that, force a rerun by hand with
  ## tar_invalidate(bench_setup_rds) once the fitting code actually
  ## changes. tests_r is still tracked/referenced below (so
  ## tar_outdated() correctly flags it when phyloslopes_tests.R changes),
  ## even though cue = "never" means that alone won't trigger a rebuild.
  tar_target(tests_r, "phyloslopes_tests.R", format = "file"),
  tar_target(
    bench_setup_rds,
    { utils_r; tests_r; source("phyloslopes_tests.R"); "phyloslopes_bench_setup.rds" },
    format = "file",
    cue = tar_cue(mode = "never")
  ),
  ## EXPENSIVE (~1000-fit microbenchmark, ~10+ min): same cue = "never"
  ## treatment as above (same caveat re: bench1_r not auto-triggering a
  ## rebuild). Reads phyloslopes_bench_setup.rds; writes two files
  ## (phylo_bench.rds / phylo_bench_slopes.rds).
  tar_target(bench1_r, "phyloslopes_bench1.R", format = "file"),
  tar_target(
    bench_rds,
    {
      utils_r; bench1_r; bench_setup_rds
      source("phyloslopes_bench1.R")
      c("phylo_bench.rds", "phylo_bench_slopes.rds")
    },
    format = "file",
    cue = tar_cue(mode = "never")
  ),
  ## fast (just plotting, no refitting) -- reruns normally when bench_rds
  ## or phyloslopes_bench1_plot.R itself changes
  tar_target(bench1_plot_r, "phyloslopes_bench1_plot.R", format = "file"),
  tar_target(
    bench_png,
    {
      bench1_plot_r; bench_rds
      source("phyloslopes_bench1_plot.R")
      c("phyloslopes_bench.png", "phyloslopes_bench1.png")
    },
    format = "file"
  ),

  ## -- write-up render --------------------------------------------------
  tar_target(qmd_source, "phyloslopes.qmd", format = "file"),
  ## phyloslopes.qmd re-derives its own data/model setup independently in
  ## its own chunks (doesn't load phyloslopes_bench_setup.rds or any tiny-
  ## example .rda), so its only real *file* dependencies are the two
  ## knitr::include_graphics() calls it makes: phyloslopes_bench1.png
  ## (tracked via bench_png, forced below) and phyloslopes_testinv_plot.png
  ## -- which is NOT wired into this pipeline (phyloslopes_testinv.R is a
  ## multi-hour, 10-core job meant to be launched by hand -- see README.md
  ## and phyloslopes_testinv.R's own header) and so must already exist on
  ## disk before this target is run. The R `quarto` package isn't
  ## installed here, so this shells out to the quarto CLI directly rather
  ## than using tarchetypes::tar_quarto().
  tar_target(
    qmd_html,
    {
      qmd_source; utils_r; bench_png
      status <- system2("quarto", c("render", "phyloslopes.qmd"))
      if (status != 0) stop("quarto render failed (status ", status, ")")
      c("phyloslopes.html", "phyloslopes.pdf")
    },
    format = "file"
  )
)
