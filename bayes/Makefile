DOCDIR=../docs/bayes

all: $(DOCDIR)/README.md $(DOCDIR)/intro.html $(DOCDIR)/bayes.bib $(DOCDIR)/priors.html $(DOCDIR)/examples.html $(DOCDIR)/setup.html

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

## also depend on examples1.rda (beyond default)
$(DOCDIR)/examples.html: examples.rmd examples1.rda
	Rscript  -e "rmarkdown::render('examples.rmd')"
	mv examples.html $(DOCDIR)

examples1.rda: run_examples.R
	R CMD BATCH --vanilla run_examples.R
	cp examples1.rda $(DOCDIR)

%.html: %.rmd
	Rscript  -e "rmarkdown::render('$<')"

test.Rout:
	Rscript -e "sink('test2.Rout'); cat(cmdstanr:::get_cmdstan_flags('STANCFLAGS'), '\n')"

## moving examples.html to docs/ not working ...

clean:
	rm -f *~ \#* *.Rout
