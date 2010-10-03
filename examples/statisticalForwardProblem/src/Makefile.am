## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/libs/core/inc \
	-I$(top_srcdir)/libs/misc/inc \
	-I$(top_srcdir)/libs/basic/inc \
	-I$(top_srcdir)/libs/stats/inc \
	-I$(TRILINOS_INCLUDE) \
	$(EPETRA_INCLUDES) \
	$(BOOST_CPPFLAGS) \
	$(GSL_CFLAGS) \
	$(HPCT_CFLAGS) \
	$(GLPK_CFLAGS)

AM_LDFLAGS  = $(BOOST_LDFLAGS) ${HPCT_LIBS} -lhpct ${GLPK_LIBS} -lglpk

LIBS        = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

noinst_PROGRAMS = exStatisticalForwardProblem_gsl

exStatisticalForwardProblem_gsl_SOURCES = exStatisticalForwardProblem_gsl.C \
	exStatisticalForwardProblem_appl.h \
	exStatisticalForwardProblem_qoi.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@