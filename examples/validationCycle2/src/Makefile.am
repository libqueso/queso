## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I. \
	      -I$(top_srcdir)/src/core/inc  \
	      -I$(top_srcdir)/src/misc/inc  \
	      -I$(top_srcdir)/src/basic/inc \
	      -I$(top_srcdir)/src/stats/inc \
	      -I$(TRILINOS_INCLUDE)          \
	        $(EPETRA_INCLUDES)           \
	        $(BOOST_CPPFLAGS)            \
	        $(GSL_CFLAGS)                \
	        $(GRVY_CFLAGS)               \
	        $(GLPK_CFLAGS)               \
	        $(HDF5_CFLAGS)

LIBS        = -L$(top_builddir)/src/ -lqueso 


exampledir  	 = @prefix@/examples/advanced
example_DATA     = tga2_gsl.C \
		   tga2_appl.C \
		   tga2_appl.h \
		   tga2_likelihood.C \
		   tga2_likelihood.h \
		   tga2_qoi.C \
		   tga2_qoi.h \
		   tga2_func.C \
		   tga2_func.h
example_PROGRAMS = tga2_gsl
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