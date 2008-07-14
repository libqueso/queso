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

OBJS_convert = black_box.eps \
	black_state.eps

all_texs: $(OBJS_convert)
	latex users.tex
	latex users.tex
	bibtex users
	latex users.tex
	latex users.tex
	dvips -t letter users.dvi -o users.ps
	ps2pdf users.ps users.pdf

.SUFFIXES: .eps .png

%.eps: %.png
	convert $< $(basename $<).eps