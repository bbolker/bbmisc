---
title: "To do"
bibliography: "bayes.bib"
---

### Examples

- Find good test for 'information gained by data' (i.e. prior/post comparison metric), e.g. https://www.princeton.edu/~umueller/priorderiv_slides.pdf: "prior sensitivity", "prior informativeness" @reimherrPrior2021; @mullerMeasuring2012

From [here](https://statmodeling.stat.columbia.edu/2019/08/10/for-each-parameter-or-other-qoi-compare-the-posterior-sd-to-the-prior-sd-if-the-posterior-sd-for-any-parameter-or-qoi-is-more-than-0-1-times-the-prior-sd-then-print-out-a-note-the-prior-dist/#comments)

> Andrew, please, please don’t put this into a publication anywhere. Even the capitulation to put it here as a reference, without qualification, isn’t a good idea. I see where you admit it’s untested but that’s the problem.<br><br>
> I’m not arguing that the idea is bad in itself. However, Cohen’s effect sizes, .05, and others, have all become unthinking “rules”, supported by references, as the thing to do. There’s very little support here for the idea and, of course, due to lack of exploration of it, no qualification of it. A better blog post for people to cite might be the exploration of this as an idea, not to focus on the 0.1 but instead to describe how one might go about arguing a prior is informative for a particular dataset in general with, at the very most, 0.1 as an example.

### Generic machinery

- how to cache results to avoid recompilation? (Can we switch priors without recompiling?)

understanding the `b_prior5`/`b_prior6` weirdness

choose priors for RE and associated FE together ...
   narrow together if regularizing
   try more intermediates, i.e. tinkering with different aspects of the prior
   
mention shinystan; Rstanarm priors
