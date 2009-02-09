## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/libs/basic/inc $(MPI_CXXFLAGS)

lib_LIBRARIES = 

lib_LIBRARIES += libuq_basic_any.a
libuq_basic_any_a_SOURCES =  uqEnvironment.C \
	uqDefines.C \
	uqVector.C \
	uqMatrix.C \
	uqMiscellaneous.C \
	uqRealFft.C \
	uqComplexFft.C \
	uqChainStatisticalOptions.C

if UQBT_GSL
lib_LIBRARIES += libuq_basic_gsl.a
libuq_basic_gsl_a_SOURCES = uqGslVector.C \
	uqGslMatrix.C \
	uqGslOdeSolver.C \
	uqGslVectorSpace.C
endif
