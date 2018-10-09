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


