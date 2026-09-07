# TODO

- [ ] Add a per-parameter error rate metric: E[V]/k per replicate (V = number of
      the k original predictors that are significant, k = original full-model
      predictor count), as an additional outcome alongside the existing
      "any significant in the model" (experimentwise) metric. For the unselected
      scenario this is trivially ~alpha (each of the k tests is independently
      exact), but for the selected scenarios it should be inflated above alpha:
      step()'s AIC-based retention and the final refit's significance check
      draw on overlapping evidence about the same predictor, so "retained" and
      "significant" are positively correlated rather than independent. Note
      this is *not* the same as literal Benjamini-Hochberg-style FDR (V/R,
      R = number of discoveries) -- since every discovery is false under this
      simulation's complete null, V/R collapses to 1 whenever R>0, making FDR
      numerically identical to the existing experimentwise (P(>=1 significant))
      metric regardless of scenario or correction method. A genuinely different
      V/R-based FDR would need true non-null effects mixed into the simulation.
- [ ] Separately, consider adding Benjamini-Hochberg and/or Benjamini-Yekutieli
      as additional *correction methods* in `screen_criteria()`, parallel to
      Holm/Dunn-Sidak/max-|T| (same "any significant after correction" outcome,
      just a different, weaker guarantee). BY is the closer analogue to Holm
      (valid under arbitrary dependence); BH is the more commonly used default
      but formally needs independence/PRDS. This is orthogonal to the
      per-parameter metric above.
