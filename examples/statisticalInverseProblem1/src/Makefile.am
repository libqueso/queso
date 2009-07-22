## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =  -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc \
	       -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(HPCT_CFLAGS)
AM_LDFLAGS  =  $(BOOST_LDFLAGS) $(BOOST_LDFLAGS) ${HPCT_LIBS} -lhpct

LIBS        = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

noinst_PROGRAMS = exStatisticalInverseProblem1_gsl
exStatisticalInverseProblem1_gsl_SOURCES = exStatisticalInverseProblem1_gsl.C exStatisticalInverseProblem1_appl.h \
	                                   exStatisticalInverseProblem1_likelihood.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@