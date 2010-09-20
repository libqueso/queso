AM_CPPFLAGS  = -I. \
	-I$(top_srcdir)/src/core/inc \
	-I$(top_srcdir)/src/misc/inc \
	-I$(top_srcdir)/src/basic/inc \
	-I$(top_srcdir)/src/stats/inc \
	-I$(TRILINOS_INCLUDE) \
	      $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(GRVY_CFLAGS) $(GLPK_CFLAGS) $(HDF5_CFLAGS)
AM_LDFLAGS   = $(BOOST_LDFLAGS) ${GRVY_LIBS} -lgrvy ${GLPK_LIBS} -lglpk ${HDF5_LIBS} -lhdf5

LIBS         = -L$(top_builddir)/src/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

EXTRA_DIST                     = TgaValidationCycle_appl.h TgaValidationCycle_likelihood.h TgaValidationCycle_qoi.h

noinst_PROGRAMS                = TgaValidationCycle_gsl
TgaValidationCycle_gsl_SOURCES = TgaValidationCycle_gsl.C 


# Required for AX_AM_MACROS
@INC_AMINCLUDE@
