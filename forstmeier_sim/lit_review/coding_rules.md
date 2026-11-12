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
