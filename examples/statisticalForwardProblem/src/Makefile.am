## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/core/inc       \
	-I$(top_srcdir)/src/misc/inc       \
	-I$(top_srcdir)/src/basic/inc      \
	-I$(top_srcdir)/src/stats/inc      \
	-I$(top_srcdir)/src/gp/inc         \
	  $(BOOST_CPPFLAGS)                \
	  $(GSL_CFLAGS)                    \
	  $(GRVY_CFLAGS)                   \
	  $(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) 
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

LIBS   = -L$(top_builddir)/src/ -lqueso


noinst_PROGRAMS = exStatisticalForwardProblem_gsl

exStatisticalForwardProblem_gsl_SOURCES = \
	exStatisticalForwardProblem_gsl.C \
	exStatisticalForwardProblem_appl.h \
	exStatisticalForwardProblem_qoi.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
