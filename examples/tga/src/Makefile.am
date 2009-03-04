## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I. -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += uqTgaExample_gsl
uqTgaExample_gsl_SOURCES = uqTgaExample_gsl.C uqTgaAppl.h uqTgaLikelihood.h uqTgaQoi.h
uqTgaExample_gsl_LDADD = $(top_builddir)/libs/libqueso.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

#bin_PROGRAMS += uqTgaExample_tri
#uqTgaExample_tri_SOURCES = uqTgaExample_tri.C uqTgaExample.h
#uqTgaExample_tri_LDADD = $(top_builddir)/libs/stats/src/libqueso_stats_gsl.a \
#	$(top_builddir)/libs/stats/src/libqueso_stats_any.a \
#	$(top_builddir)/libs/basic/src/libqueso_basic_gsl.a \
#	$(top_builddir)/libs/basic/src/libqueso_basic_tri.a \
#	$(top_builddir)/libs/basic/src/libqueso_basic_any.a \
#	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@