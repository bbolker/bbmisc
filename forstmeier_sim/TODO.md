# TODO

- [ ] Decide whether to add `bonferroni` to the criteria (regex and
      `criterion_labels`) in `R/plot_corrections_by_scenario_core.R`
- [ ] cosmetics: github link/icon?
- [ ] Consider adding Benjamini-Hochberg and/or Benjamini-Yekutieli
      as additional *correction methods* in `screen_criteria()`, parallel to
      Holm/Dunn-Sidak/max-|T| (same "any significant after correction" outcome,
      just a different, weaker guarantee). BY is the closer analogue to Holm
      (valid under arbitrary dependence); BH is the more commonly used default
      but formally needs independence/PRDS.
- [ ] Consider switching Holm for Bonferroni in Figure 2 (J. Dushoff's
      suggestion)

