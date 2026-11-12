# TODO

- [ ] Decide whether to add `bonferroni` to the criteria (regex and
      `criterion_labels`) in `R/plot_corrections_by_scenario_core.R`

- [ ] Sort out Figure 2 (`output/fig1_corrections.png`). For the
      experimentwise error rate, Holm is identical to Bonferroni (Holm
      rejects something iff min p < alpha/k), and Dunn-Sidak is nearly so.
      Options:
    - put the uncorrected panel back in and drop Holm
    - switch to the per-parameter version
      (`output/fig1_corrections_perparam.png`), where Holm does differ
    - reduce to two panels: Bonferroni and single-step max|T|

- [ ] cosmetics: github link/icon?
- [ ] Consider adding Benjamini-Hochberg and/or Benjamini-Yekutieli
      as additional *correction methods* in `screen_criteria()`, parallel to
      Holm/Dunn-Sidak/max-|T| (same "any significant after correction" outcome,
      just a different, weaker guarantee). BY is the closer analogue to Holm
      (valid under arbitrary dependence); BH is the more commonly used default
      but formally needs independence/PRDS.

