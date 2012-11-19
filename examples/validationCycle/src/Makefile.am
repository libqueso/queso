## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I. \
	      -I$(top_srcdir)/src/core/inc  \
	      -I$(top_srcdir)/src/misc/inc  \
	      -I$(top_srcdir)/src/basic/inc \
	      -I$(top_srcdir)/src/stats/inc \
	        $(BOOST_CPPFLAGS)           \
	        $(GSL_CFLAGS)               \
	        $(GRVY_CFLAGS)              \
		$(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

LIBS        = -L$(top_builddir)/src/ -lqueso

# 11/16/2012: No longer creating examples in the installed QUESO directory
#exampledir  	 = @prefix@/examples/advanced
#example_DATA     = exTgaValidationCycle_gsl.C exTgaValidationCycle_appl.h \
#		   exTgaValidationCycle_likelihood.h exTgaValidationCycle_qoi.h
#example_PROGRAMS = exTgaValidationCycle_gsl
#
#exTgaValidationCycle_gsl_SOURCES = exTgaValidationCycle_gsl.C exTgaValidationCycle_appl.h \
#				   exTgaValidationCycle_likelihood.h exTgaValidationCycle_qoi.h

noinst_PROGRAMS = exTgaValidationCycle_gsl

exTgaValidationCycle_gsl_SOURCES = exTgaValidationCycle_gsl.C \
								exTgaValidationCycle_appl.h \
								exTgaValidationCycle_likelihood.h \
								exTgaValidationCycle_qoi.h

#bin_PROGRAMS += exTgaValidationCycle_tri
#exTgaValidationCycle_tri_SOURCES = exTgaValidationCycle_tri.C exTgaValidationCycle_appl.h exTgaValidationCycle_likelihood.h exTgaValidationCycle_qoi.h
#exTgaValidationCycle_tri_LDADD = $(top_builddir)/libs/stats/src/libqueso_stats_gsl.a \
#	$(top_builddir)/libs/stats/src/libqueso_stats_any.a \
#	$(top_builddir)/libs/basic/src/libqueso_basic_gsl.a \
#	$(top_builddir)/libs/basic/src/libqueso_basic_tri.a \
#	$(top_builddir)/libs/basic/src/libqueso_basic_any.a \
#	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
