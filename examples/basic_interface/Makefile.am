## Process this file with automake to produce Makefile.in

#exampledir		= @prefix@/examples/basic
#example_DATA            = inverse_problem.c inverse_problemf.f90 queso.inp

#
EXTRA_DIST              = inverse_problemf.f90 inverse_problem.h queso.inp

AM_CFLAGS               = $(GRVY_CFLAGS)
AM_CPPFLAGS             = -I$(top_srcdir)/src/interface
LIBS                    = -L$(top_builddir)/src/ -lqueso

#example_PROGRAMS        = inverse_problem latency

#inverse_problem_SOURCES = inverse_problem.c
#latency_SOURCES         = latency_model.c

#aqui
noinst_PROGRAMS = inverse_problem latency

inverse_problem_SOURCES = inverse_problem.c
latency_SOURCES         = latency_model.c

@INC_AMINCLUDE@

