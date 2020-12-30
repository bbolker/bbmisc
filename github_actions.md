---
title: "stuff about GitHub actions"
date: 28 December 2020
---

## introduction

Many of my R packages are set up using **continuous integration** to run tests every time changes are pushed to a repository. I have been using Travis-CI, but they are phasing out their free-tier support for open source projects, so (like much of the R community) I am migrating to GitHub Actions (GA). (I don't know what options are available for Gitlab/Bitbucket and other Git repositories that are not GitHub ...)

## general information

- At present GA is [free for public repositories](https://docs.github.com/en/free-pro-team@latest/actions/reference/usage-limits-billing-and-administration), there are limits (2000 minutes/month in the free tier) for private repos.
- GA is controlled by a YAML configuration file (**much** more about this below); this is also true of Travis builds, but I don't remember fighting with Travis configuration being quite so painful (perhaps it was and I've just repressed it/forgotten about it).
- there is of course lots of [documentation](https://docs.github.com/en/free-pro-team@latest/actions) about GA, but I'm lazy and tried to jump right in ...
- most of what I did is based on the [usethis package](https://usethis.r-lib.org/), and in particular on [Dean Attali's blog post on migrating from Travis to GA](https://deanattali.com/blog/migrating-travis-to-github). In the simplest case, you "just" need to:
   - un-enable building on Travis's website
   - remove Travis artifacts from your package (especially `.travis.yml`)
   - run the appropriate `usethis::` command (see below for details)
   - optionally, add an appropriate badge URL to your README (automatically suggested by `usethis` machinery
   - et voilà, you're done!
   
**BUT** ...

## fussy bits

- the `usethis` functionality relating to GA does not seem super well-documented. `use_github_actions()` is the simplest/default way to set everything up, but `use_github_action_check_standard()` is more thorough; it's what Dean Attali recommends in his blog post
   - `use_github_actions()` sets up a much smaller/simpler YAML configuration. It checks *only* on MacOS (why MacOS??? why not Linux?), and configures the testing platform with a smaller number of built-in system packages (especially, LaTeX support). I configured one package (`McMasterPandemic`) with this version.
       - added stuff: [pandoc](https://github.com/r-lib/actions/tree/master/setup-pandoc), tinytex (I should have used the [r-lib action](https://github.com/r-lib/actions/tree/master/setup-tinytex) but did it by hand)
   - `use_github_action_check_standard()` 
       - seems to automatically incorporate `texlive` (better/more complete than `tinytex`?)
       - but tests on 4 platforms (MacOS, Linux-release, Linux-devel, Windows) automatically, which is overkill for my needs: I used it anyway (for its better coverage) and commented out the platforms I didn't want
	   - supposedly caches packages (but not seeing much on Linux-devel?)

### working directories

If your package is in a subdirectory you need to sprinkle in appropriate `working-directory:` specifications

### skipping tests

Travis-CI has a built-in test that automatically skips testing when your commit message contains the string "[skip ci]", and lists the tests as having been skipped. There is a [way to do this in GA](https://github.com/marketplace/actions/ci-skip-action), but (while shiny and sophisticated) it doesn't do what I want, i.e. easily allow skipping/not skipping in a way that tags the skipped builds as 'skipped' rather than 'failed'. 

Constantin Ahlmann-Eltze suggested on twitter that I look at this solution from the `matrixStats` package:

```r
if: "! contains(github.event.head_commit.message, '[ci skip]')"
```

I am going to use this with `"[skip ci]"` instead (already familar from Travis-CI). This approach also work with an "opt-in" statement (for my colleague who feels that CI is contributing to planetary degradation):

```r
if: "contains(github.event.head_commit.message, '[run ci]')"
```

This appears to work as intended (in `lme4`). When tests are skipped they show up as **skipped** rather than failed. It would be nice to work out a convenient/pre-packaged way to choose among options (skip, standard, thorough/multi-platform) based on tags in the commit messages.


### YAML format

- It seems *incredibly* fussy about whitespace etc. Tabs instead of spaces will break it. A one-character error in whitespace will break it. (GA reports the line on which the error occurred). I used emacs `whitespace-mode` to help, but it would be good if there were a linter that I knew how to install (to be able to get the syntax right *before* pushing/discovering it's broken on GH).
- You need to know the difference between a line starting with a dash (the beginning of a list) and without one (an item within a list). The all-important `steps:` section alternates sections that are either `- uses:` (which runs the contents of a GH repo/page/branch somewhere, e.g. to apply a recipe for setting up `pandoc`) or `- name:` (which names a step to run). Within these steps we have **non-dash-prefixed lines at the same indentation level as the tag** such as `working-directory:` or `run:` that actually specify what to do/how to do it.

### LaTeX stuff

- building LaTeX vignettes is quite likely to fail
- need either `tinytex` (R package + install machinery) or `texlive` installed
    - extra packages for `tinytex`: `tinytex::tlmgr_install("foo")` (where `foo` is the LaTeX package name) in a step that uses `shell: Rscript {0}`
	- extra packages for `texlive`: `sudo apt-get install ...` in a regular (shell-script) step. In future I might just add `sudo apt-get install texlive-science texlive-latex-extra texlive-bibtex-extra` by default

- Spent a long time finding problematic packages in LaTeX headers one at a time. This helps if using `texlive`: 
```r
debsrch <- function(s) {
    browseURL(sprintf("https://packages.debian.org/search?suite=default&section=all&arch=any&searchon=contents&keywords=%s",s))
}
```
- If using `tinytex`, `tinytex::tlmgr_search("something.sty")` in a local R session works well to figure out what LaTeX package you need

### Miscellaneous

If you're going to experiment with this stuff on a new branch, make sure to add your branch name to the `branches:` section at the top of the file (otherwise nothing will happen)

### My packages

Hopefully these are useful examples.

- [McMasterPandemic](https://github.com/bbolker/McMasterPandemic/blob/master/.github/workflows/R-CMD-check.yaml)
   - basic (`uses_github_actions()`)
   - added `"[skip ci]"` à la `matrixStats`
   - used `setup-pandoc` from `r-lib/actions`
   - added `tinytex` by hand, plus a whole pile of `tinytex::tlmgr_install()` statements
       - could have combined into a single long character vector?
	   - could install `texlive` + some extra Debian packages instead?
- [glmmTMB](https://github.com/glmmTMB/glmmTMB/blob/github_actions/.github/workflows/R-CMD-check.yaml)
   - used `uses_github_actions_standard()`
   - commented out extra platforms, left just `R-release`
   - added texlive packages
   - modified working directories (code gratefully accepted from Bryce Mecum https://github.com/amoeba)
- [lme4](https://github.com/lme4/lme4/blob/master/.github/workflows/R-CMD-check.yaml)
   - used `..._standard()`, commented out all but `r-devel` platform
   - added texlive: `sudo apt-get install texlive texlive-latex-base texlive-latex-extra`
   - added `build_args="--compact-vignettes=both"` to the `rcmdcheck` command ([no quotation marks around "both"](https://stat.ethz.ch/pipermail/r-package-devel/2020q4/006099.html))
	   - installing dependencies is very slow on `r-devel` before the first successful test, but after that package caching appears to work as intended
   - opt-in to tests via `"[run ci]"` in commit messages
- [fitode](https:://github.com/parksw3/fitode/blog/master/.github/workflows/R-CMD-check.yaml)
   - similar to those above: `standard()` with other platforms commented, `[skip ci]`
   - no need for subdirectories, extra vignette compacting, etc.
   - added a step to install `texlive` and friends, but not sure I really needed to
   - in the end I found that a custom `Makefile` in the `vignettes/` directory was messing things up (generating a "vignette product does not have a known filename extension" error), even though it was OK on my machine *and* was included in `.Rbuildignore`?? 
**to do**: `bbmle`, ?
