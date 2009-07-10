all:	all_texs

clean: clean_texs clean_exec

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

OBJS_convert = figs/black_box.eps \
	figs/black_state.eps \
	figs/black_box3.eps \
	figs/himmelAobs.eps \
	figs/himmelAcomp.eps \
	figs/himmel_compMH_chain1.eps \
	figs/himmel_compMH_chain1_2d.eps \
	figs/himmel_compDR_chain1.eps \
	figs/himmel_compDR_chain1_2d.eps \
	figs/himmel_compAM_chain1.eps \
	figs/himmel_compAM_chain1_2d.eps \
	figs/himmel_compDRAM_chain1.eps \
	figs/himmel_compDRAM_chain1_2d.eps

#all_texs: $(OBJS_convert)
all_texs: 
	latex users.tex
	latex users.tex
	bibtex users
	latex users.tex
	latex users.tex
	dvips -t letter users.dvi -o users.ps
	ps2pdf users.ps users_0_4_0.pdf

.SUFFIXES: .eps .png

%.eps: %.png
	convert $< $(basename $<).eps