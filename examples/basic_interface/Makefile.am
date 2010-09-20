## Process this file with automake to produce Makefile.in

exampledir		= @prefix@/examples/basic
example_DATA            = inverse_problem.c inverse_problemf.f90 queso.inp

EXTRA_DIST              = inverse_problemf.f90 inverse_problem.h queso.inp

AM_CFLAGS               = $(GRVY_CFLAGS)
AM_CPPFLAGS             = -I$(top_srcdir)/src/interface
AM_LDFLAGS              = $(BOOST_LDFLAGS) ${GRVY_LIBS} -lgrvy ${GLPK_LIBS} -lglpk ${HDF5_LIBS} -lhdf5
LIBS                    = -L$(top_builddir)/src/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS) ${HDF5_LIBS} -lhdf5

example_PROGRAMS        = inverse_problem latency

inverse_problem_SOURCES = inverse_problem.c
latency_SOURCES         = latency_model.c

@INC_AMINCLUDE@
