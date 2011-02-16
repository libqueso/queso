## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/core/inc       \
	-I$(top_srcdir)/src/misc/inc       \
	-I$(top_srcdir)/src/basic/inc      \
	-I$(top_srcdir)/src/stats/inc      \
	  $(BOOST_CPPFLAGS)                \
	  $(GSL_CFLAGS)                    \
	  $(GRVY_CFLAGS)                   \
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

LIBS   = -L$(top_builddir)/src/ -lqueso $(ANN_LIBS)

noinst_PROGRAMS = exInfoTheory_gsl

exInfoTheory_gsl_SOURCES = \
	exInfoTheory_gsl.C \
	exInfoTheory_appl.h
#	exInfoTheory_qoi.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@