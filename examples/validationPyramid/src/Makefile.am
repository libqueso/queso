## Process this file with automake to produce Makefile.in

DIST_SUBDIRS = physics_1 physics_2 physics_1_2

AM_CPPFLAGS = -I. \
	      -I$(top_srcdir)/examples/validationPyramid/src/physics_1   \
	      -I$(top_srcdir)/examples/validationPyramid/src/physics_2   \
	      -I$(top_srcdir)/examples/validationPyramid/src/physics_1_2 \
	      -I$(top_srcdir)/src/core/inc  \
	      -I$(top_srcdir)/src/misc/inc  \
	      -I$(top_srcdir)/src/basic/inc \
	      -I$(top_srcdir)/src/stats/inc \
	      -I$(TRILINOS_INCLUDE) \
	        $(EPETRA_INCLUDES)  \
	        $(BOOST_CPPFLAGS)   \
	        $(GSL_CFLAGS)       \
	        $(GRVY_CFLAGS)      \
	        $(GLPK_CFLAGS)      \
	        $(HDF5_CFLAGS)

LIBS        = -L$(top_builddir)/src/ -lqueso 

noinst_PROGRAMS                 = exValidationPyramid_gsl
exValidationPyramid_gsl_SOURCES = exValidationPyramid_gsl.C


# Required for AX_AM_MACROS
@INC_AMINCLUDE@
