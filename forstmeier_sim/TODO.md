# TODO

- [ ] Consider reporting error rates restricted to just main-effect terms
      (excluding interactions), separate from the existing all-terms
      metrics -- e.g. a `n_sig_raw_main_*`-style count over only the m
      terms whose name has no `:`, alongside the existing count over all k
      terms. Motivation: now that main effects are tested at the
      population mean via sum-to-zero contrasts, it's worth checking
      whether they carry the same selection-bias inflation as interaction
      terms under `step()`, or a different degree of it (`step()`'s
      marginality constraint means an interaction's survival depends on
      its main effects' presence but not vice versa, so the two term types
      needn't behave symmetrically). Corrected criteria (Dunn-Sidak/Holm/
      max-|T|) would presumably still calibrate on the full k (the actual
      multiple-testing exposure) but only check significance among the
      main-effect subset. The current output only stores the all-terms
      aggregate (`n_sig_raw_*`, not a per-term breakdown by main-effect vs.
      interaction), so this would need a `screen_criteria()`/
      `one_rep_corrections()` schema extension and another full simulation
      rerun, the same pattern as the n_sig_* addition.
- [ ] Consider adding Benjamini-Hochberg and/or Benjamini-Yekutieli
      as additional *correction methods* in `screen_criteria()`, parallel to
      Holm/Dunn-Sidak/max-|T| (same "any significant after correction" outcome,
      just a different, weaker guarantee). BY is the closer analogue to Holm
      (valid under arbitrary dependence); BH is the more commonly used default
      but formally needs independence/PRDS.
