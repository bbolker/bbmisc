DOCDIR=../docs/bayes

all: $(DOCDIR)/README.md $(DOCDIR)/intro.html $(DOCDIR)/bayes.bib $(DOCDIR)/priors.html $(DOCDIR)/examples.html

## bizarrely, examples.rmd fails when run via 'make'
##   but NOT when run via Rscript -e ... from the command line???
## Error in if (stancflags_local != "") { : the condition has length > 1
## Calls: <Anonymous> ... eval -> .fun -> <Anonymous> -> initialize -> <Anonymous>

## fix: one specific rule for moving html to $DOCDIR, another more general rule for copying to $DOCDIR
$(DOCDIR)/bayes.bib: bayes.bib
	mkdir -p $(DOCDIR)
	cp $< $(DOCDIR)

$(DOCDIR)/README.md: README.md
	mkdir -p $(DOCDIR)
	cp $< $(DOCDIR)


$(DOCDIR)/%: % bayes.bib
	mkdir -p $(DOCDIR)
	mv $< $(DOCDIR)

%.html: %.rmd
	Rscript  -e "rmarkdown::render('$<')"


