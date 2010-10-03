## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/libs/core/inc \
	-I$(top_srcdir)/libs/misc/inc       \
	-I$(top_srcdir)/libs/basic/inc      \
	-I$(top_srcdir)/libs/stats/inc      \
	-I$(TRILINOS_INCLUDE)               \
	  $(EPETRA_INCLUDES)                \
	  $(BOOST_CPPFLAGS)                 \
	  $(GSL_CFLAGS)                     \
	  $(GRVY_CFLAGS)                    \
	  $(GLPK_CFLAGS)                    \
	  $(HDF5_CFLAGS)

AM_LDFLAGS  = $(BOOST_LDFLAGS) ${GRVY_LIBS} -lqueso ${GLPK_LIBS} -lglpk ${HDF5_LIBS} -lhdf5
LIBS        = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

noinst_PROGRAMS = exStatisticalForwardProblem_gsl

exStatisticalForwardProblem_gsl_SOURCES = exStatisticalForwardProblem_gsl.C \
	exStatisticalForwardProblem_appl.h \
	exStatisticalForwardProblem_qoi.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@