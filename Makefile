## This is bbmisc

current: target
-include target.mk

######################################################################

## Content

rmdfiles = $(wildcard *.rmd)
rmdh = $(rmdfiles:%.rmd=%.html)

Sources += $(rmdfiles)
Ignore += $(rmdh)

peak_I.html: peak_I.rmd
peak_I_simple.html: peak_I_simple.rmd

peak_reduction.Rout: peak_reduction.R

######################################################################

## Spline stuff 2022 Nov 07 (Mon)

## git diff dca3b2c0 -- Rmisc/spline_quantiles.R
Sources += Rmisc/spline_quantiles.R
Rmisc/spline_quantiles.Rout: Rmisc/spline_quantiles.R
	$(pipeR)

######################################################################

## Rules from Bolker

%.html: %.qmd
	quarto render $<

docs/%.html: %.html
	mv $< docs/$<

%.pdf: %.pdf
	quarto render $< --to pdf

%.html: %.[Rr]md
	Rscript -e "library(\"rmarkdown\"); render(\"$<\")"

## https://bbolker.github.io/bbmisc/brant_survive.html
## brant_survive.html: brant_survive.rmd

%.html: %.md
	Rscript -e "library(\"rmarkdown\"); render(\"$<\")"

%.md: %.rmd
	Rscript -e "library(\"knitr\"); knit(\"$<\")"

%.tex: %.Rnw
	Rscript -e "library(\"knitr\"); knit(\"$<\")"

%.tex: %.md
	pandoc -s -S -t latex -V documentclass=tufte-handout $*.md -o $*.tex

%.pdf: %.md
	echo "rmarkdown::render(\"$<\", output_format=\"pdf_document\")" | R --slave

%.pdf: %.tex
	pdflatex --interaction=nonstopmode $*

clean:
	rm -f *.log *.aux *.md *.out *.nav *.snm *.toc *.vrb texput.log *~

peeves: peeves.md
	Rscript -e 'rmarkdown::render("peeves.md")'
	mv peeves.html docs/

rtips: r_parallel_hpc.rmd
	Rscript -e 'rmarkdown::render("r_parallel_hpc.rmd")'
	mv r_parallel_hpc.html docs/

######################################################################

## makestuff needs to be made manually, since we want to allow Bolker-style making without makestuff
alldirs += bayes

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff
makestuff: makestuff/Makefile
makestuff/Makefile:
	(ls ../makestuff/Makefile && /bin/ln -s ../makestuff) || git clone $(msrepo)/makestuff
	ls $@

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
