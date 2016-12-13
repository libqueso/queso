## Process this file with automake to produce Makefile.in

DIST_SUBDIRS = physics_1 physics_2 physics_1_2

AM_CPPFLAGS =
AM_CPPFLAGS += -I.
AM_CPPFLAGS += -I$(top_srcdir)/examples/validationPyramid/src/physics_1
AM_CPPFLAGS += -I$(top_srcdir)/examples/validationPyramid/src/physics_2
AM_CPPFLAGS += -I$(top_srcdir)/examples/validationPyramid/src/physics_1_2
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)
AM_CPPFLAGS += $(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) 
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

LIBS        = -L$(top_builddir)/src/ -lqueso 

noinst_PROGRAMS                 = exValidationPyramid_gsl
exValidationPyramid_gsl_SOURCES = exValidationPyramid_gsl.C


# Required for AX_AM_MACROS
@INC_AMINCLUDE@
