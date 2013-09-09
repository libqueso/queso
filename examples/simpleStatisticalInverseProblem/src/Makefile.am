## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += -I$(top_builddir)/src/core/inc  # For queso.h
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)
AM_CPPFLAGS += $(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) 
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

LIBS        = -L$(top_builddir)/src/ -lqueso $(HDF5_LIBS)

noinst_PROGRAMS = exSimpleStatisticalInverseProblem_gsl

exStatisticalInverseProblem_gsl_SOURCES = \
	example_main.C \
	example_compute.C \
	example_likelihood.C \
	example_compute.h \
	example_likelihood.h
	       
# Required for AX_AM_MACROS
@INC_AMINCLUDE@
