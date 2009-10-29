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

OBJS_convert = figs/uqEnvironment.eps \
	figs/uqEnvironmentOptions.eps \
	figs/uqVectorSet.eps \
	figs/uqVector.eps \
	figs/uqMatrix.eps \
	figs/uqScalarFunction.eps \
	figs/uqVectorFunction.eps \
	figs/uqScalarSequence.eps \
	figs/uqVectorSequence.eps \
	figs/uqSip.eps \
	figs/uqSipOptions.eps \
	figs/uqSfp.eps \
	figs/uqSfpOptions.eps \
	figs/uqMetropolisHastingsSG.eps \
	figs/uqMetropolisHastingsSGOptions.eps \
	figs/uqMonteCarloSG.eps \
	figs/uqMonteCarloSGOptions.eps \
	figs/uqSequenceStatisticalOptions.eps \
	figs/uqVectorRandomVariable.eps \
	figs/paper_plot1.eps \
	figs/paper_plot2.eps \
	figs/paper_plot3.eps \
	figs/paper_plot4.eps \
	figs/paper_plot5.eps \
	figs/queso_paper1_03.eps \
	figs/queso_paper1_09.eps \
	figs/queso_paper1_10.eps

all_texs: $(OBJS_convert)
	latex users.tex
	latex users.tex
	bibtex users
	latex users.tex
	latex users.tex
	dvips -t letter users.dvi -o users.ps
	ps2pdf users.ps users_0_41_0.pdf

.SUFFIXES: .eps .png .dia

%.eps: %.png
	convert $< $(basename $<).eps

%.eps: %.dia
	dia -t eps $<
