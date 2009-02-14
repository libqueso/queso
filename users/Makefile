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
	black_state.eps \
	black_box3.eps \
	himmelAobs.eps \
	himmelAcomp.eps \
	himmel_compMH_chain1.eps \
	himmel_compMH_chain1_2d.eps \
	himmel_compDR_chain1.eps \
	himmel_compDR_chain1_2d.eps \
	himmel_compAM_chain1.eps \
	himmel_compAM_chain1_2d.eps \
	himmel_compDRAM_chain1.eps \
	himmel_compDRAM_chain1_2d.eps

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