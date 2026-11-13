# Coding rules

Decisions on what counts as "a GLM with two or more explanatory variables"
(F&S 2011), applied when filling in `coding.csv`.

- Mixed models (LMM, GLMM) count.
- GAMs and Bayesian regressions do not count; record them in `model_class`.
- RNA-seq and other genomics analyses (e.g. per-gene differential-expression
  GLMs in DESeq2/edgeR/limma, GWAS) do not count. A paper whose only models are
  of this kind is coded `uses_glm_2plus = n` and is dropped from the
  denominator of empirical studies; if it also has ordinary GLMs, code those.
- Phylogenetic GLS (PGLS) and other frequentist regressions with correlated
  errors count as LMs, provided they have two or more predictors.

## Full model presented

- `y` if results for every term of the full (global, pre-selection) model are
  reported anywhere, including the supplement: coefficient tables, or a full
  ANOVA/Wald/LRT table of every term (even without parameter estimates).
- A paper that never simplifies and reports its fitted model counts as `y`.
- An AIC/AICc model-comparison table that lists the full model as one of the
  candidates does not count, since it gives no estimates for its terms.

## Model selection

- Choosing the random-effects structure by AIC (or dropping random effects by
  any criterion) does not count as model selection; only data-driven choice
  of fixed-effect terms (predictors or interactions) does. Code such papers
  `model_selection = none (...)`, with the random-effects step described in
  parentheses.
  (@matuschekBalancingTypeError2017 recommend backward selection of random
  effects with a relaxed p-value criterion, and find that it barely inflates
  type I error.)

## Multiple-comparison correction

- `mcc` records corrections for the number of model terms (or models) tested,
  e.g. an experiment-wide Holm correction.
- Post hoc pairwise adjustments within a factor (Tukey, emmeans-adjusted
  contrasts, Bonferroni/Holm on pairwise comparisons) go in
  `posthoc_adjustment` instead; they do not count as `mcc`.
- In either column, a correction applied to a model after selection counts
  only as "minimal k": it corrects for the retained terms, not for those in
  the original model. Add "minimal k: applied after selection" inside the
  parentheses.

## Numeric summaries

- `n_explanatory_max`, `n_terms_max`, `N_over_k_min`: one number per paper,
  from the largest model (predictors; terms including interactions) and the
  model with the smallest N/k. `Inf` = "large"; blank = unknown; "<2" is
  stored as 2. The free-text columns keep the details.

## Treatment x covariate interactions

- `trt_x_covariate` (experimental or "both" designs only): an interaction
  between a manipulated treatment and a variable that was measured rather
  than assigned (sex, population, species, site, date). Treatment x treatment
  interactions do not count. `y?` marks borderline cases, which are counted
  as `y`.
