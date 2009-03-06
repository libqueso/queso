## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I. -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += exTgaValidationCycle_gsl
exTgaValidationCycle_gsl_SOURCES = exTgaValidationCycle_gsl.C exTgaValidationCycle_appl.h exTgaValidationCycle_likelihood.h exTgaValidationCycle_qoi.h
exTgaValidationCycle_gsl_LDADD = $(top_builddir)/libs/libqueso.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

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