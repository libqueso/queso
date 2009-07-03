## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/libs/interface
AM_LDFLAGS  = $(BOOST_LDFLAGS) -L${HPCT_LIB} -lhpct

bin_PROGRAMS            = inverse_problem
inverse_problem_SOURCES = inverse_problem.c
inverse_problem_LDADD   = $(top_builddir)/libs/libqueso.a  \
		          $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)
@INC_AMINCLUDE@