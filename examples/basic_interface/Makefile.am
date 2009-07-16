## Process this file with automake to produce Makefile.in

exampledir		= @prefix@/examples
example_DATA            = inverse_problem.c queso.inp
#EXTRA_DIST              = inverse_problem.c queso.inp

AM_CPPFLAGS             = -I$(top_srcdir)/libs/interface
AM_LDFLAGS              = $(BOOST_LDFLAGS) ${HPCT_LIBS} -lhpct
LIBS                    = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

example_PROGRAMS        = inverse_problem
inverse_problem_SOURCES = inverse_problem.c

@INC_AMINCLUDE@