# Wald vs LRT: what's known, and reconciling the results

Most of the individual pieces are known, but they're spread across several literatures that rarely cite each other. What's less common is seeing them in one place. The belief that "the LRT is always preferable" is folk wisdom that goes further than the theory supports. All the references below are from memory and haven't been checked against the full texts. Several aren't in the document's bib file, so they're worth checking before relying on them.

## What's well established

**1. The Wald statistic fails for large effects (Hauck-Donner).** This is the classic result: Hauck & Donner (1977), Væth (1985), Fears et al. (1996), Pawitan (2000), Yee (2022). The Wald statistic isn't monotone in the effect size, so it loses power and gives absurd intervals near separation. The references in the document cover this well.

**2. Small-sample intervals for the odds ratio are hard, and discreteness makes coverage oscillate.**

- **Binomial intervals:** Brown, Cai & DasGupta (2001, *Statistical Science*) describe the erratic, sawtooth coverage of binomial intervals, likelihood-ratio ones included.
- **Odds ratios:** Agresti (1999, *Biometrics*, "On logit confidence intervals for the odds ratio with small samples") shows that adding pseudo-counts improves the logit Wald interval. That's the Haldane–Anscombe/Gart adjustment, which is the Firth Wald interval in the document.
- **Systematic comparisons:** the closest match to the binomial section is probably Fagerland, Lydersen & Laake (2015, *Statistical Methods in Medical Research*), "Recommended confidence intervals for two independent binomial proportions", along with their book *Statistical Analysis of Contingency Tables* (2017). As I recall, they compare many odds-ratio intervals in this exact two-group setting. Their recommendations are neither plain Wald nor plain profile likelihood, but intervals such as Baptista-Pike mid-p and the adjusted inverse sinh.

**3. Firth/Jeffreys penalization: finite estimates, and shrinkage at large effects.**

- **Firth (1993):** introduces the method, which removes the O(1/n) bias.
- **Heinze & Schemper (2002):** apply it to separation, using penalized profile-likelihood intervals.
- **Puhr et al. (2017, *Statistics in Medicine*):** show that Firth biases predicted probabilities towards ½ with rare events.
- **Kosmidis & Firth (2021, *Biometrika*), "Jeffreys-prior penalty, finiteness and shrinkage in binomial-response GLMs":** as I recall, this is the formal statement of the "overcorrection at large β". The Jeffreys penalty shrinks estimates towards zero.
- **Greenland & Mansournia (2015, *Statistics in Medicine*):** discuss log-F priors as pseudo-data alternatives to Firth. That's the same kind of construction as the t-prior approximation.

**4. The theory that resolves most of the puzzle.** The main theoretical case for the LRT is:

- **Invariance:** it gives the same answer however the parameter is written; Wald does not.
- **Following the likelihood's shape:** it isn't forced to be symmetric.
- **Better two-sided accuracy:** two-sided LR intervals have coverage error O(1/n) and can be Bartlett-corrected.

But the **signed root** of the LR statistic, which governs each tail separately, has error O(1/√n), the same order as Wald. The two tail errors largely cancel in the two-sided interval. This is standard in higher-order asymptotics. It's the motivation for the modified signed root r\* (Barndorff-Nielsen; Pierce & Peters 1992; Brazzale, Davison & Reid 2007, *Applied Asymptotics*). So the finding that the LRT is fine two-sided but unbalanced between the tails is predicted by theory.

## Reconciling this with "LRT is preferable"

Where the LRT reliably wins is **power and sensible behaviour of tests**, and **interval shape** under strong asymmetry: no Hauck-Donner collapse, and no symmetric intervals crossing impossible values. Nothing guarantees that its **coverage** always beats Wald's, and in two places in the results it doesn't:

- **Moderate β (about 1–4):** here the unpenalized Wald looks better, because two errors cancel. The MLE's upward bias pushes the interval up, and the too-narrow symmetric shape pulls it back. The LRT keeps the bias but fixes the shape, so its lower-tail errors show up. The penalized comparisons confirm this: once the bias is removed, the profile interval beats the Wald interval at every β. So the comparison says LRT beats Wald when the estimator is well behaved. The plain MLE isn't well behaved here, and Wald hides that by accident.
- **Large β:** Wald "wins" on coverage by giving up and reporting (−∞, ∞) under separation. That's coverage bought with zero information, which is the flip side of the Hauck-Donner power loss. A method with 100% coverage and no power isn't better.

## How unreliable is the LRT, really?

Less than the plots may suggest, because the failures are concentrated in a narrow region:

- **Where things are typical** (β up to about 3 at n ≥ 50, i.e. odds ratios up to about 20 with baseline probability 0.2): two-sided coverage is within about 0.01–0.02 of nominal for every method.
- **The unpenalized LRT's worst two-sided coverage** is about 0.91. It occurs in the band where separation has become common but the true β is still below the fixed lower limit that separation produces: β ≈ 5 for n = 100, β ≈ 3.5 for n = 20. With baseline 0.2, β = 5 means p₁ ≈ 0.97, which is an extreme setting.
- **The penalized profile intervals** keep two-sided coverage at about 0.93–0.97 for n ≥ 50. Their remaining problem is the tail imbalance at large β, from overcorrection.

## What's less standard in the document

These parts are, as far as I know, less commonly written up:

- **Splitting coverage by tail and relating each tail to bias vs interval shape:** the higher-order-asymptotics literature makes this point in the abstract, but seldom with pictures like these.
- **Measuring separation's contribution to the LRT's undercoverage:** most of the lower-tail miss at β ≈ 4 for n = 50 comes from separated outcomes.
- **A t prior available in glmmTMB that closely approximates Firth, placed on the group logits.**
- **Penalized profile vs penalized Wald comparisons that separate bias from interval shape.**

Papers most worth downloading, to check these claims and whether the results are already covered:

- Fagerland et al. (2015)
- Agresti (1999)
- Kosmidis & Firth (2021)
- Heinze & Schemper (2002)
- Brazzale, Davison & Reid (2007) or Pierce & Peters (1992), for the signed-root point

It may also be worth softening the "LRT preferable" framing in the tl;dr along these lines. Something like: "the LRT's advantages are in power and interval shape; its one-sided coverage is not more accurate than Wald's in general, and when the estimator is biased the Wald interval can look better by accident."
