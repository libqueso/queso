all:	all_texs

clean_all: clean_texs clean_exec

clean_texs:
	rm -f *.aux
	rm -f *.log
	rm -f *.dvi
	rm -f *.toc
	rm -f *.blg
	rm -f *.bbl
	rm -f users*.ps
	rm -f *.pdf

clean_exec:
	rm -f *~

all_texs:
	latex users.tex
	latex users.tex
	bibtex users
	latex users.tex
	latex users.tex
	dvips -t letter users.dvi -o users.ps
	ps2pdf users.ps users.pdf
