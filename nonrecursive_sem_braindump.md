# non-recursive structural equation modeling

What do I (BMB) know about NRSEM?

## SEM in general

This is a method for partitioning effects occuring along different causal pathways. Shipley (2002) is a readable introduction. Originally from Sewall Wright's method of *path analysis*. At its simplest, you can do a path analysis by (1) writing down a path diagram, i.e. hypothesizing a particular set of causal pathways; (2) running a series of linear models that quantify proportions of variation that act along different pathways. In this way you can (in principle) separate variation explained according to different direct and indirect pathways.

* everything is easier if you assume that all relationships are linear and all of observed variables are Gaussian
* you _may_ be able to test the validity of different models
* results depend _strongly_ on your willingness to commit to a specific model (e.g. ruling out the effects of other exogenous variables; if exogenous variables are affecting everything in the system you're probably screwed)
* important to distinguish _exploratory_ from _confirmatory_ SEM.
* [lavaan](https://CRAN.R-project.org/package=lavaan) is probably the most popular R package, see also [sem](https://cran.r-project.org/package=sem), less friendly but more powerful (has 2SLS for NRSEM ...). [Mplus](https://www.statmodel.com/) (commercial) is probably the most powerful tool available (there are R interfaces).
* measurement error makes everything harder

## nonrecursive models

* recursive models have *uncorrelated* error terms and *unidirectional* causation (no loops or reciprocal paths); non-recursive models have correlated error terms and/or loops/reciprocal paths
* *block recursive* models may have loops etc. within specific blocks/subsets of variables
* *identifiability* concerns: can we actually estimate the parameters?
* *instrumental variables* are extra observations that may fix the problem
* need fancier methods, e.g. *two-stage least squares* or *full-information maximum likelihood* (both implemented in the `sem` package)
* could work through the Duncan, Haller and Portes example: see `?sem::effects.sem` and the Fox Barcelona notes

## worrying about causality

**but** the most important questions are probably "why am I doing this? Do I trust that my model is complete and correctly specified? Does it pass the 'smell test' that I would plausibly be able to derive these kinds of conclusions from this kind of data?"

Doing this with *cross-sectional data* (individuals measured once) is very, very, tricky. Estimating causality is hard even when we have longitudinal/time series data (see [Granger causality](https://en.wikipedia.org/wiki/Granger_causality)) ("does not capture instantaneous and non-linear causal relationships")

Also Sugihara et al 2012; Cobey and Baskerville 2016 ("Although CCM has theoretical support, natural systems routinely violate its assumptions"); Sugihara et al 2017 (response to C&B 2016); Correia et al 2025

## miscellaneous

Connections of DAGs (directed acyclic graphs) to analysis of moderators, confounders, colliders, etc., see McElreath *Rethinking Statistics*.

## Links

* https://stats.stackexchange.com/questions/577602/sem-recursive-versus-nonrecursive-models

>  Non-recursive models are very rare and (in my experience) they have all kinds of problems - both theoretical (Really? You think that's meaningful with cross-sectional data) and practical (they don't converge). IMHO they should just have a footnote that says "This is sometimes possible, but probably don't do it."

* [Fox Barcelona notes](https://facsocsci.mcmaster.ca/jfox/Courses/R/IQSBarcelona/SEMs-notes.pdf)
* [CRAN Psychometrics Task View](https://cran.r-project.org/web/views/Psychometrics.html)
* https://en.wikipedia.org/wiki/Structural_equation_modeling

 Shipley, Bill. 2002. Cause and Correlation in Biology: A User’s Guide to Path Analysis, Structural Equations and Causal Inference. 1st ed. Cambridge University Press.

Cobey, Sarah, and Edward B. Baskerville. 2016. “Limits to Causal Inference with State-Space Reconstruction for Infectious Disease.” PLOS ONE 11 (12): e0169050. https://doi.org/10.1371/journal.pone.0169050.

Correia, Hannah E., Laura E. Dee, Jarrett E. K. Byrnes, et al. 2025. Best Practices for Moving from Correlation to Causation in Ecological Research. June 6. https://ecoevorxiv.org/repository/view/9361/.

Sugihara, George, Ethan R. Deyle, and Hao Ye. 2017. “Reply to Baskerville and Cobey: Misconceptions about Causation with Synchrony and Seasonal Drivers.” Proceedings of the National Academy of Sciences 114 (12): E2272–74. https://doi.org/10.1073/pnas.1700998114.

Sugihara, George, Robert May, Hao Ye, et al. 2012. “Detecting Causality in Complex Ecosystems.” Science (New York, N.Y.) 338 (6106): 496–500. https://doi.org/10.1126/science.1227079.
