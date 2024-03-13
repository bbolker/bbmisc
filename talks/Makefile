%.html: %.rmd
	echo "rmarkdown::render(\"$<\")" | R --slave

%.html: %.qmd
	quarto render $<

