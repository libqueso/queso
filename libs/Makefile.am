## Process this file with automake to produce Makefile.in

# Build in these directories:
# DIST_SUBDIRS = basic stats

AM_CPPFLAGS = -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc \
	      -I$(top_srcdir)/libs/interface $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(HPCT_CFLAGS)


# Dynamic library build

lib_LTLIBRARIES     = libqueso.la
libqueso_la_LDFLAGS = $(all_libraries) -release $(GENERIC_RELEASE)

# Sources from basic/src

libqueso_la_SOURCES = $(top_srcdir)/libs/basic/src/uqEnvironment.C \
	$(top_srcdir)/libs/basic/src/uqEnvironmentOptions.C \
	$(top_srcdir)/libs/basic/src/uqDefines.C \
	$(top_srcdir)/libs/basic/src/uqVector.C \
	$(top_srcdir)/libs/basic/src/uqMatrix.C \
	$(top_srcdir)/libs/basic/src/uqMiscellaneous.C \
	$(top_srcdir)/libs/basic/src/uqRealFft.C \
	$(top_srcdir)/libs/basic/src/uqComplexFft.C \
	$(top_srcdir)/libs/basic/src/uqSequenceStatisticalOptions.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVector.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosMatrix.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVectorSpace.C



# Sources from basic/src added conditional on gsl
if UQBT_GSL

libqueso_la_SOURCES += $(top_srcdir)/libs/basic/src/uqGslVector.C       \
 	               $(top_srcdir)/libs/basic/src/uqGslMatrix.C       \
 	               $(top_srcdir)/libs/basic/src/uqGslVectorSpace.C  \
 	               $(top_srcdir)/libs/basic/src/uqGslAsciiTable.C
endif

# Sources for Basic Interface 

libqueso_la_SOURCES += $(top_srcdir)/libs/interface/basic_classes.cpp   \
                      $(top_srcdir)/libs/interface/basic_interface.cpp

#libqueso_a_LIBADD   = $(HPCT_LIBSTATIC)

# Sources from stats/src

libqueso_la_SOURCES += $(top_srcdir)/libs/stats/src/uqMarkovChainSGOptions.C \
	$(top_srcdir)/libs/stats/src/uqMLSamplingOptions.C \
	$(top_srcdir)/libs/stats/src/uqMLSamplingLevelOptions.C \
	$(top_srcdir)/libs/stats/src/uqMonteCarloSGOptions.C \
	$(top_srcdir)/libs/stats/src/uqStatisticalInverseProblemOptions.C \
	$(top_srcdir)/libs/stats/src/uqStatisticalForwardProblemOptions.C

libqueso_includedir = $(prefix)/include

# Headers to install from basic/inc
libqueso_include_HEADERS = $(top_srcdir)/libs/basic/inc/uq1D1DFunction.h \
	$(top_srcdir)/libs/basic/inc/uq2dArrayOfStuff.h \
	$(top_srcdir)/libs/basic/inc/uqArrayOfOneDGrids.h \
	$(top_srcdir)/libs/basic/inc/uqArrayOfOneDTables.h \
	$(top_srcdir)/libs/basic/inc/uqArrayOfSequences.h \
	$(top_srcdir)/libs/basic/inc/uqAsciiTable.h \
	$(top_srcdir)/libs/basic/inc/uqSequenceStatisticalOptions.h \
	$(top_srcdir)/libs/basic/inc/uqDefines.h \
	$(top_srcdir)/libs/basic/inc/uqEnvironment.h \
	$(top_srcdir)/libs/basic/inc/uqEnvironmentOptions.h \
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
	$(top_srcdir)/libs/stats/inc/uqMarkovChainSGOptions.h \
	$(top_srcdir)/libs/stats/inc/uqMLSampling.h \
	$(top_srcdir)/libs/stats/inc/uqMLSamplingOptions.h \
	$(top_srcdir)/libs/stats/inc/uqMLSamplingLevelOptions.h \
	$(top_srcdir)/libs/stats/inc/uqModelValidation.h \
	$(top_srcdir)/libs/stats/inc/uqMonteCarloSG.h \
	$(top_srcdir)/libs/stats/inc/uqMonteCarloSGOptions.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalForwardProblem.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalForwardProblemOptions.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalInverseProblem.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalInverseProblemOptions.h \
	$(top_srcdir)/libs/stats/inc/uqScalarCdf.h \
	$(top_srcdir)/libs/stats/inc/uqTKGroup.h \
	$(top_srcdir)/libs/stats/inc/uqValidationCycle.h \
	$(top_srcdir)/libs/stats/inc/uqVectorCdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorMdf.h \
	$(top_srcdir)/libs/stats/inc/uqJointPdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorRealizer.h \
	$(top_srcdir)/libs/stats/inc/uqVectorRV.h

# Headers to install from Basic Interface

libqueso_include_HEADERS += $(top_srcdir)/libs/interface/queso.h

