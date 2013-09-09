## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += -I$(top_builddir)/src/core/inc  # For queso.h
AM_CPPFLAGS += -I$(TRILINOS_INCLUDE)
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)
AM_CPPFLAGS += $(GLPK_CFLAGS)
AM_CPPFLAGS += $(HDF5_CFLAGS)

LIBS        = -L$(top_builddir)/src/ -lqueso 

noinst_PROGRAMS = exStatisticalInverseProblem1_gsl

exStatisticalInverseProblem1_gsl_SOURCES = \
	exStatisticalInverseProblem1_gsl.C \
	exStatisticalInverseProblem1_appl.h \
	exStatisticalInverseProblem1_likelihood.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
