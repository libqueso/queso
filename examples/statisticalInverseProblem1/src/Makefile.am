## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =  -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc \
	 -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += uqStatisticalInverseProblem1_gsl
uqStatisticalInverseProblem1_gsl_SOURCES = uqStatisticalInverseProblem1_gsl.C uqStatisticalInverseProblem1_appl.h uqStatisticalInverseProblem1_likelihood.h
uqStatisticalInverseProblem1_gsl_LDADD = $(top_builddir)/libs/libqueso.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@