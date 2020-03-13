
## This is bbmisc

current: target
-include target.mk


######################################################################

%.html: %.rmd
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

######################################################################

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff
makestuff: makestuff/Makefile
makestuff/Makefile:
	(ls ../makestuff/Makefile && /bin/ln -s ../makestuff) || git clone $(msrepo)/makestuff
	ls $@

-include makestuff/os.mk

## -include makestuff/wrapR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

