# Lit review TODO

## 2026 survey (all 45 sampled papers coded; supplements fetched where needed)

- [x] Settle definitions of full model / model selection / MCC (now in
  `coding_rules.md`) and recode.
- [ ] Spot-check the remaining medium-confidence codings in `coding.csv`
  against the papers, applying the definitions in `coding_rules.md`.
- [ ] Mention random-effects selection in the main discussion
  (`forstmeier_sig.qmd`), citing Matuschek et al. 2017.
- [x] Summary table: numeric `*_max`/`*_min` columns, treatment x covariate
  row, exact binomial CIs.
- [x] Genomics-only papers dropped from the denominator (not replaced).
- [ ] Check the two borderline `trt_x_covariate = y?` codes (ecy.70478,
  rspb.2026.0555).

## 2007 calibration (after 2026 is done)

- [ ] Rerun the same pipeline (TOC, abstract screening, sampling, fetching,
  coding with `coding_rules.md`) on the September 2007 issues of the same five
  journals and compare with F&S (28/50 = 56%, 6/28 with 6-17 predictors, 80%
  interactions, 16/16 simplified, 2/28 full models). Keep outputs separate
  (`*_sept2007.csv`). Decide whether to code all empirical papers or sample
  10 per journal.

## Notes

- Downloaded papers live in `papers/` (gitignored); a copy is in
  `~/Dropbox/forstmeier_sim/papers`.
- Animal Behaviour (Elsevier) PDFs must be downloaded by hand; Elsevier
  supplements are on the public CDN
  (`ars.els-cdn.com/content/image/1-s2.0-<PII>-mmc1.<ext>`).
