## Process this file with automake to produce Makefile.in

DIST_SUBDIRS = physics_1 physics_2 physics_1_2

AM_CPPFLAGS = -I. -I$(top_srcdir)/examples/pyramid/src/physics_1 \
	-I$(top_srcdir)/examples/pyramid/src/physics_2 \
	-I$(top_srcdir)/examples/pyramid/src/physics_1_2 \
	-I$(top_srcdir)/libs/basic/inc \
	-I$(top_srcdir)/libs/stats/inc \
	-I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += uqPyramidExample_gsl
uqPyramidExample_gsl_SOURCES = uqPyramidExample_gsl.C
uqPyramidExample_gsl_LDADD = $(top_builddir)/libs/libqueso.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
