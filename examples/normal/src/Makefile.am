## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =  -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc \
	 $(MPI_CXXFLAGS) $(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
AM_LDFLAGS = $(BOOST_LDFLAGS)

bin_PROGRAMS = 

bin_PROGRAMS += uqNormalExample_gsl
uqNormalExample_gsl_SOURCES = uqNormalExample_gsl.C uqNormalAppl.h uqNormalLikelihood.h
uqNormalExample_gsl_LDADD = $(top_builddir)/libs/libqueso.a \
	$(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(MPI_CXXLIBS) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@