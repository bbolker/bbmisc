%.html: %.rmd
	echo "rmarkdown::render('$*.rmd')" | R --slave

%.pdf: %.Rnw
	echo "knitr::knit2pdf('$*.Rnw')" | R --slave
