DOCDIR=../docs/bayes

all: $(DOCDIR)/README.md $(DOCDIR)/intro.html $(DOCDIR)/bayes.bib

## fix: one specific rule for moving html to $DOCDIR, another more general rule for copying to $DOCDIR
$(DOCDIR)/bayes.bib: bayes.bib
	mkdir -p $(DOCDIR)
	cp $< $(DOCDIR)

$(DOCDIR)/README.md: README.md
	mkdir -p $(DOCDIR)
	cp $< $(DOCDIR)


$(DOCDIR)/%: %
	mkdir -p $(DOCDIR)
	mv $< $(DOCDIR)

%.html: %.rmd
	Rscript  -e "rmarkdown::render('$<')"


