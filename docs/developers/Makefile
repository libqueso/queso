all:	all_texs

clean_all: clean_texs clean_exec

clean_texs:
	rm -f *.aux
	rm -f *.log
	rm -f *.dvi
	rm -f *.toc
	rm -f *.blg
	rm -f *.bbl
	rm -f dev*.ps
	rm -f *.pdf

clean_exec:
	rm -f *~

all_texs:
	latex dev.tex
	latex dev.tex
	bibtex dev
	latex dev.tex
	latex dev.tex
	dvips -t letter dev.dvi -o dev.ps
	ps2pdf dev.ps developers.pdf
