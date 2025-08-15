## Specific goals

* evaluate the frequency of use of two specific different analytical methods over time/journal/etc.. Specifically, we're considering the use of generalized linear models (GLMs) with either *power-logistic* or *complementary log-log* links to evaluate the effects of covariates on the probability of nest predation in birds
* develop/illustrate the use of open-source/local tools, and automation, to do these kinds of analyses
* see [here](https://bbolker.github.io/bbmisc/brant_survive.html) for some background, and Shaffer papers

## Paper selection

* Literature search with appropriate keywords: "avian nest predation" is pretty good
* This can be done with commercial tools (Web of Science [available via Mac libraries], Scopus, ...) but I would like to try OpenAlex.
* e.g. https://openalex.org/works?page=1&filter=title_and_abstract.search:avian+nest+predation gets 924 hits, which seems like a good starting point; 780 are articles
* (I've found OpenAlex's metadata dodgy in the past, but it might not hurt for this example)
* adding 'power-logistic' *in title/abstract* gets only 4 ... https://openalex.org/works?page=1&filter=title_and_abstract.search:avian+nest+predation+power-logistic

## Paper downloads

* automate downloads of PDFs, conversions to text
* not sure how tricky it will be to download papers that are not open source but are licensed to Mac (authentication/cookies/etc.)

## Paper screening

* want to subdivide papers according to which method they use (or neither).
* need some ground-truth/human evaluations
* can we use ChatGPT (see Mitchell and Earn paper)
* can we use a local, open-source LLM such as Ollama?
* need to develop prompts (see Earn and Mitchell)

## tool development

* can we/should we do 'vibe coding', i.e. use LLMs to build the code? (At least we have a clear idea of what we want the machinery to do ...)

## references

Mitchell, Evan, Elisha B. Are, Caroline Colijn, and David J. D. Earn. 2025. “Using Artificial Intelligence Tools to Automate Data Extraction for Living Evidence Syntheses.” PLOS ONE 20 (4): e0320151. https://doi.org/10.1371/journal.pone.0320151.

Keck, Francois, Henry Broadbent, and Florian Altermatt. 2025. “Extracting Massive Ecological Data on State and Interactions of Species Using Large Language Models.” Preprint, bioRxiv, January 27. https://doi.org/10.1101/2025.01.24.634685.

Moorthy, Sruthi M. Krishna, Man Qi, Alice Rosen, Yadvinder Malhi, and Rob Salguero-Gomez. 2025. Harnessing Large Language Models for Ecological Literature Reviews: A Practical Pipeline. February 7. https://ecoevorxiv.org/repository/view/8516/.

Gougherty, Andrew V., and Hannah L. Clipp. 2024. “Testing the Reliability of an AI-Based Large Language Model to Extract Ecological Information from the Scientific Literature.” Npj Biodiversity 3 (1): 13. https://doi.org/10.1038/s44185-024-00043-9.


Rotella, Jay J., Stephen J. Dinsmore, and Terry L. Shaffer. 2004. “Modeling Nest-Survival Data: A Comparison of Recently Developed Methods That Can Be Implemented in MARK and SAS.” Animal Biodiversity and Conservation 27 (1): 187–205.

Shaffer, Terry L. 2004. “A Unified Approach to Analyzing Nest Success.” The Auk 121 (2): 526–40. https://doi.org/10.1093/auk/121.2.526.

Heisey, DENNIS M., TERRY L. Shaffer, and GARY C. White. 2007. “The ABCs of Nest Survival: Theory and Application from a Biostatistical Perspective.” Studies in Avian Biology 34: 13.
