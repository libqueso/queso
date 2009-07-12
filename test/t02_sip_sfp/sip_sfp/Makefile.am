AM_CPPFLAGS  = -I. -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc -I$(TRILINOS_INCLUDE) \
	      $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
AM_LDFLAGS   = $(BOOST_LDFLAGS)

EXTRA_DIST                = example_compute.h example_likelihood.h example_qoi.h

bin_PROGRAMS              = SipSfpExample_gsl
SipSfpExample_gsl_SOURCES = example_main.C example_compute.C example_likelihood.C example_qoi.C
SipSfpExample_gsl_LDADD   = $(top_builddir)/libs/libqueso.a \
				 $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
