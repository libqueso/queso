## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I. -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/queso/inc $(MPI_CXXFLAGS) $(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += uqTgaExample
uqTgaExample_SOURCES = uqTgaExample_gsl.C uqTgaExample.h
uqTgaExample_LDADD = $(top_builddir)/libs/queso/src/libuq_queso_gsl.a \
	$(top_builddir)/libs/queso/src/libuq_queso_any.a \
	$(top_builddir)/libs/basic/src/libuq_basic_gsl.a \
	$(top_builddir)/libs/basic/src/libuq_basic_any.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(MPI_CXXLIBS) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@