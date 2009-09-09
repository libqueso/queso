AM_CPPFLAGS  = -I. \
	-I$(top_srcdir)/libs/core/inc \
	-I$(top_srcdir)/libs/misc/inc \
	-I$(top_srcdir)/libs/basic/inc \
	-I$(top_srcdir)/libs/stats/inc \
	-I$(TRILINOS_INCLUDE) \
	      $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(HPCT_CFLAGS)
AM_LDFLAGS   = $(BOOST_LDFLAGS) ${HPCT_LIBS} -lhpct

EXTRA_DIST   = example_compute.h example_likelihood.h example_qoi.h
LIBS         = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

noinst_PROGRAMS           = SipSfpExample_gsl
SipSfpExample_gsl_SOURCES = example_main.C example_compute.C example_likelihood.C example_qoi.C

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
