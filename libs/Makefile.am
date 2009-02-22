## Process this file with automake to produce Makefile.in

# Build in these directories:
# DIST_SUBDIRS = basic stats

AM_CPPFLAGS = -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc $(MPI_CXXFLAGS) $(BOOST_CPPFLAGS)

lib_LIBRARIES = libqueso.a

# Sources from basic/src
libqueso_a_SOURCES = $(top_srcdir)/libs/basic/src/uqEnvironment.C \
	$(top_srcdir)/libs/basic/src/uqDefines.C \
	$(top_srcdir)/libs/basic/src/uqVector.C \
	$(top_srcdir)/libs/basic/src/uqMatrix.C \
	$(top_srcdir)/libs/basic/src/uqMiscellaneous.C \
	$(top_srcdir)/libs/basic/src/uqRealFft.C \
	$(top_srcdir)/libs/basic/src/uqComplexFft.C \
	$(top_srcdir)/libs/basic/src/uqChainStatisticalOptions.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVector.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosMatrix.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVectorSpace.C

# Sources from basic/src added conditional on gsl
if UQBT_GSL
libqueso_a_SOURCES += $(top_srcdir)/libs/basic/src/uqGslVector.C \
	$(top_srcdir)/libs/basic/src/uqGslMatrix.C \
	$(top_srcdir)/libs/basic/src/uqGslVectorSpace.C
endif

# Sources from stats/src
# Currently none... add if necessary

libqueso_includedir = $(prefix)/include

# Headers to install from basic/inc
libqueso_include_HEADERS = $(top_srcdir)/libs/basic/inc/uq1D1DFunction.h \
	$(top_srcdir)/libs/basic/inc/uq2dArrayOfStuff.h \
	$(top_srcdir)/libs/basic/inc/uqArrayOfOneDGrids.h \
	$(top_srcdir)/libs/basic/inc/uqArrayOfOneDTables.h \
	$(top_srcdir)/libs/basic/inc/uqArrayOfSequences.h \
	$(top_srcdir)/libs/basic/inc/uqAsciiTable.h \
	$(top_srcdir)/libs/basic/inc/uqChainStatisticalOptions.h \
	$(top_srcdir)/libs/basic/inc/uqDefines.h \
	$(top_srcdir)/libs/basic/inc/uqEnvironment.h \
	$(top_srcdir)/libs/basic/inc/uqFft.h \
	$(top_srcdir)/libs/basic/inc/uqGslMatrix.h \
	$(top_srcdir)/libs/basic/inc/uqGslVector.h \
	$(top_srcdir)/libs/basic/inc/uqInstantiateIntersection.h \
	$(top_srcdir)/libs/basic/inc/uqMatrix.h \
	$(top_srcdir)/libs/basic/inc/uqMiscellaneous.h \
	$(top_srcdir)/libs/basic/inc/uqOneDGrid.h \
	$(top_srcdir)/libs/basic/inc/uqScalarFunction.h \
	$(top_srcdir)/libs/basic/inc/uqScalarFunctionSynchronizer.h \
	$(top_srcdir)/libs/basic/inc/uqScalarSequence.h \
	$(top_srcdir)/libs/basic/inc/uqSequenceOfVectors.h \
	$(top_srcdir)/libs/basic/inc/uqTrilinosMatrix.h \
	$(top_srcdir)/libs/basic/inc/uqTrilinosVector.h \
	$(top_srcdir)/libs/basic/inc/uqVectorFunction.h \
	$(top_srcdir)/libs/basic/inc/uqVectorFunctionSynchronizer.h \
	$(top_srcdir)/libs/basic/inc/uqVector.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSequence.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSet.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSpace.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSubset.h

# Headers to install from stats/inc
libqueso_include_HEADERS += $(top_srcdir)/libs/stats/inc/uqCovCond.h \
	$(top_srcdir)/libs/stats/inc/uqMarkovChainPositionData.h \
	$(top_srcdir)/libs/stats/inc/uqMarkovChainSG1.h \
	$(top_srcdir)/libs/stats/inc/uqMarkovChainSG2.h \
	$(top_srcdir)/libs/stats/inc/uqModelValidation.h \
	$(top_srcdir)/libs/stats/inc/uqMonteCarloSG.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalForwardProblem.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalInverseProblem.h \
	$(top_srcdir)/libs/stats/inc/uqScalarCdf.h \
	$(top_srcdir)/libs/stats/inc/uqTKGroup.h \
	$(top_srcdir)/libs/stats/inc/uqValidationCycle.h \
	$(top_srcdir)/libs/stats/inc/uqVectorCdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorMdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorPdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorRealizer.h \
	$(top_srcdir)/libs/stats/inc/uqVectorRV.h

