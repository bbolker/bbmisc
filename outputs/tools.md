# development tools we should be using more


## auto-formatters/linters

I still format inconsistently (part of this is based on some emacs
configuration junk I haven’t bothered to sort out). We can use *linters*
to detect format inconsistencies (and other kinds of statically
detectable code problems), but it would be/could be more efficient to
use *auto-formatters* (along with whatever formatting aids are built
into your editor/IDE) to clean this stuff up and minimize the number of
distracting code changes that are only from different
formatting/whitespace/etc. by different contributors.
[air](https://posit-dev.github.io/air/formatter.html) is a fast
formatter (from Posit); unfortunately, it’s *opinionated*: \> Air is
purposefully minimally configurable, with the main configuration points
being related to line width and indent style. … so, for example, it
won’t allow leading commas (see
https://bsky.app/profile/usrbinr.bsky.social/post/3lz2s24gzd22q ). Maybe
we could vibe-code something for `styler`, see
[here](https://styler.r-lib.org/reference/create_style_guide.html) ?

## environment preservation

There are a bunch of tools for this: `renv`,
[uvr](https://nbafrank.github.io/uvr/), Docker … like the
auto-formatters, part of the problem is the degree of
opinionatedness/how much these tools impose themselves on all users of a
system (I’ve had some hassles with `renv` …)

## 

Quarto, typst

## AI tools

Of course.

## workflow tools

make, shellpipes/makestuff, targets

## IDEs

VSCode, Positron

## testing frameworks

lazytest, tinytest, snapshots

## 

Are we relying too heavily on GitHub? Proprietary, more and more
AI-embedded. While the core git functionality makes sure that we can
migrate to another platform/git server any time we want, how do we make
sure that the ‘arounds’ (issues, pull request discussions, etc.) aren’t
locked in?

https://pub.towardsai.net/finding-a-way-out-a-deep-dive-into-github-alternatives-in-2026-3c13de3226e8

Codeberg? https://github.com/git-bug/git-bug ? What about actions/CI-CD?
