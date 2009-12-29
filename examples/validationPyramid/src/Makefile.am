## Process this file with automake to produce Makefile.in

DIST_SUBDIRS = physics_1 physics_2 physics_1_2

AM_CPPFLAGS = -I. \
	-I$(top_srcdir)/examples/validationPyramid/src/physics_1 \
	-I$(top_srcdir)/examples/validationPyramid/src/physics_2 \
	-I$(top_srcdir)/examples/validationPyramid/src/physics_1_2 \
	-I$(top_srcdir)/libs/core/inc \
	-I$(top_srcdir)/libs/misc/inc \
	-I$(top_srcdir)/libs/basic/inc \
	-I$(top_srcdir)/libs/stats/inc \
	-I$(TRILINOS_INCLUDE) \
	$(EPETRA_INCLUDES) \
	$(BOOST_CPPFLAGS) \
	$(GSL_CFLAGS) \
	$(HPCT_CFLAGS) \
	$(GLPK_CFLAGS)

AM_LDFLAGS  = $(BOOST_LDFLAGS) ${HPCT_LIBS} -lhpct ${GLPK_LIBS} -lglpk

LIBS        = -L$(top_builddir)/libs/.libs -lqueso $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(EPETRA_LIBS)

noinst_PROGRAMS                 = exValidationPyramid_gsl
exValidationPyramid_gsl_SOURCES = exValidationPyramid_gsl.C


# Required for AX_AM_MACROS
@INC_AMINCLUDE@
