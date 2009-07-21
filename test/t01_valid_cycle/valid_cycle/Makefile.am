AM_CPPFLAGS  = -I. -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc -I$(TRILINOS_INCLUDE) \
	      $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(HPCT_CFLAGS)
AM_LDFLAGS   = $(BOOST_LDFLAGS) ${HPCT_LIBS} -lhpct

LIBS         = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

EXTRA_DIST                     = TgaValidationCycle_appl.h TgaValidationCycle_likelihood.h TgaValidationCycle_qoi.h

bin_PROGRAMS                   = TgaValidationCycle_gsl
TgaValidationCycle_gsl_SOURCES = TgaValidationCycle_gsl.C 


# Required for AX_AM_MACROS
@INC_AMINCLUDE@
