## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src/core/inc  \
	      -I$(top_srcdir)/src/misc/inc  \
	      -I$(top_srcdir)/src/basic/inc \
	      -I$(top_srcdir)/src/stats/inc \
	      -I$(top_srcdir)/src/gp/inc    \
	      -I$(TRILINOS_INCLUDE) \
	         \
	        $(BOOST_CPPFLAGS) \
	        $(GSL_CFLAGS) \
	        $(GRVY_CFLAGS) \
	        $(GLPK_CFLAGS) \
	        $(HDF5_CFLAGS)

LIBS        = -L$(top_builddir)/src/ -lqueso 

noinst_PROGRAMS = exStatisticalInverseProblem1_gsl

exStatisticalInverseProblem1_gsl_SOURCES = \
	exStatisticalInverseProblem1_gsl.C \
	exStatisticalInverseProblem1_appl.h \
	exStatisticalInverseProblem1_likelihood.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
