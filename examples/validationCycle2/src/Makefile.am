## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I. \
	      -I$(top_srcdir)/src/core/inc  \
	      -I$(top_srcdir)/src/misc/inc  \
	      -I$(top_srcdir)/src/basic/inc \
	      -I$(top_srcdir)/src/stats/inc \
	        $(BOOST_CPPFLAGS)            \
	        $(GSL_CFLAGS)                \
	        $(GRVY_CFLAGS)               \
		$(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

LIBS        = -L$(top_builddir)/src/ -lqueso 

# 11/16/2012: No longer creating examples in the installed QUESO directory
#exampledir  	 = @prefix@/examples/advanced
#example_DATA     = tga2_gsl.C \
#		   tga2_appl.C \
#		   tga2_appl.h \
#		   tga2_likelihood.C \
#		   tga2_likelihood.h \
#		   tga2_qoi.C \
#		   tga2_qoi.h \
#		   tga2_func.C \
#		   tga2_func.h
#example_PROGRAMS = tga2_gsl
#tga2_gsl_SOURCES = tga2_gsl.C \
#		   tga2_appl.C \
#		   tga2_appl.h \
#		   tga2_likelihood.C \
#		   tga2_likelihood.h \
#		   tga2_qoi.C \
#		   tga2_qoi.h \
#		   tga2_func.C \
#		   tga2_func.h

noinst_PROGRAMS = tga2_gsl
tga2_gsl_SOURCES = tga2_gsl.C \
		   tga2_appl.C \
		   tga2_appl.h \
		   tga2_likelihood.C \
		   tga2_likelihood.h \
		   tga2_qoi.C \
		   tga2_qoi.h \
		   tga2_func.C \
		   tga2_func.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
