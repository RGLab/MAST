# Makefile to use knitr for package vignettes

# put all PDF targets here, separated by spaces
PDFS=SingleCellAssay-intro.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *.dvi *tikzDictionary 

%.pdf: %.vignette
	cp symbols.foo symbols.tex
	cp $*.vignette $*.Rnw
	$(R_HOME)/bin/R  -e "library(knitr); knit2pdf('$*.Rnw')"
	cp $*.dummy $*.Rnw
	rm symbols.tex
