## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src/core/inc \
	-I$(top_srcdir)/src/misc/inc       \
	-I$(top_srcdir)/src/basic/inc      \
	-I$(top_srcdir)/src/stats/inc      \
	-I$(TRILINOS_INCLUDE)               \
	  $(EPETRA_INCLUDES)                \
	  $(BOOST_CPPFLAGS)                 \
	  $(GSL_CFLAGS)                     \
	  $(GRVY_CFLAGS)                    \
	  $(GLPK_CFLAGS)                    \
	  $(HDF5_CFLAGS)

AM_LDFLAGS  = $(BOOST_LDFLAGS) ${GRVY_LIBS} -lqueso ${GLPK_LIBS} -lglpk ${HDF5_LIBS} -lhdf5
LIBS        = -L$(top_builddir)/src/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

noinst_PROGRAMS = exStatisticalForwardProblem_gsl

exStatisticalForwardProblem_gsl_SOURCES = exStatisticalForwardProblem_gsl.C \
	exStatisticalForwardProblem_appl.h \
	exStatisticalForwardProblem_qoi.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@