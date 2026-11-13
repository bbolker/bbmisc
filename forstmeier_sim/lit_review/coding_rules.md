# Coding rules

Decisions on what counts as "a GLM with two or more explanatory variables"
(F&S 2011), applied when filling in `coding.csv`.

- Mixed models (LMM, GLMM) count.
- GAMs and Bayesian regressions do not count; record them in `model_class`.
- RNA-seq and other genomics analyses (e.g. per-gene differential-expression
  GLMs in DESeq2/edgeR/limma, GWAS) do not count. A paper whose only models are
  of this kind is coded `uses_glm_2plus = n` but stays in the denominator of
  empirical studies; if it also has ordinary GLMs, code those.
- Phylogenetic GLS (PGLS) and other frequentist regressions with correlated
  errors count as LMs, provided they have two or more predictors.

## Full model presented

- Coefficient tables for the full (global, pre-selection) model count,
  including when they appear only in the supplement.
- A full ANOVA table (Type II/III, or Wald chi-squared tests) for every term in
  the full model also counts, even without parameter estimates.

## Model selection

- Choosing the random-effects structure by AIC (or dropping random effects by
  any criterion) does not count as model selection; only selection among
  fixed-effect terms does. Code such papers `model_selection = none (...)`,
  with the random-effects step described in parentheses.
  (@matuschekBalancingTypeError2017 recommend backward selection of random
  effects with a relaxed p-value criterion, and find that it barely inflates
  type I error.)

## Multiple-comparison correction (mcc)

- Post hoc adjustments (Tukey, emmeans-adjusted contrasts, Bonferroni/Holm on
  pairwise comparisons) count.
- A correction applied to a model after selection counts only as "minimal k":
  it corrects for the retained terms, not for those in the original model.
  Add "minimal k: applied after selection" inside the `mcc` parentheses.
- An AIC/AICc model-comparison table that lists the full model as one of the
  candidates does not count, since it gives no estimates for its terms.
