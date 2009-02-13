## Process this file with automake to produce Makefile.in

DIST_SUBDIRS = physics_1 physics_2 physics_1_2

AM_CPPFLAGS = -I. -I$(top_srcdir)/examples/queso/pyramid/src/physics_1 \
	-I$(top_srcdir)/examples/queso/pyramid/src/physics_2 \
	-I$(top_srcdir)/examples/queso/pyramid/src/physics_1_2 \
	-I$(top_srcdir)/libs/basic/inc \
	-I$(top_srcdir)/libs/stats/inc \
	$(MPI_CXXFLAGS) $(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += uqPyramidExample
uqPyramidExample_SOURCES = uqPyramidExample_gsl.C
uqPyramidExample_LDADD = $(top_builddir)/libs/libqueso.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(MPI_CXXLIBS) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
