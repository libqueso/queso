## Process this file with automake to produce Makefile.in

# Build in these directories:
# DIST_SUBDIRS = basic stats

AM_CPPFLAGS = -I$(top_srcdir)/libs/core/inc \
	-I$(top_srcdir)/libs/misc/inc \
	-I$(top_srcdir)/libs/basic/inc \
	-I$(top_srcdir)/libs/stats/inc \
	-I$(top_srcdir)/libs/interface \
	$(BOOST_CPPFLAGS) \
	$(GSL_CFLAGS) \
	$(HPCT_CFLAGS) \
	$(GLPK_CFLAGS)


# Dynamic library build

lib_LTLIBRARIES     = libqueso.la
libqueso_la_LDFLAGS = $(all_libraries) -release $(GENERIC_RELEASE)

# (1/5) Sources from core/src
libqueso_la_SOURCES = $(top_srcdir)/libs/core/src/uqDefines.C \
	$(top_srcdir)/libs/core/src/uqEnvironment.C \
	$(top_srcdir)/libs/core/src/uqEnvironmentOptions.C \
	$(top_srcdir)/libs/core/src/uqVector.C \
	$(top_srcdir)/libs/core/src/uqMatrix.C \
	$(top_srcdir)/libs/core/src/uqMPI.C \
	$(top_srcdir)/libs/core/src/uqTrilinosVector.C \
	$(top_srcdir)/libs/core/src/uqTrilinosMatrix.C

# (1/5) Sources from core/src added conditional on gsl
if UQBT_GSL
libqueso_la_SOURCES += $(top_srcdir)/libs/core/src/uqGslVector.C \
	$(top_srcdir)/libs/core/src/uqGslMatrix.C
endif

# (2/5) Sources from misc/src
libqueso_la_SOURCES += $(top_srcdir)/libs/misc/src/uq1D1DFunction.C \
	$(top_srcdir)/libs/misc/src/uq1DQuadrature.C \
	$(top_srcdir)/libs/misc/src/uqComplexFft.C \
	$(top_srcdir)/libs/misc/src/uqMiscellaneous.C \
	$(top_srcdir)/libs/misc/src/uqRealFft.C

# (2/5) Sources from misc/src added conditional on gsl
if UQBT_GSL
libqueso_la_SOURCES += $(top_srcdir)/libs/misc/src/uqGslAsciiTable.C
endif

# (3/5) Sources from basic/src
libqueso_la_SOURCES += $(top_srcdir)/libs/basic/src/uqSequenceStatisticalOptions.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVectorSpace.C

# (3/5) Sources from basic/src added conditional on gsl
if UQBT_GSL
libqueso_la_SOURCES += $(top_srcdir)/libs/basic/src/uqGslVectorSpace.C
endif

# (4/5) Sources from stats/src
libqueso_la_SOURCES += $(top_srcdir)/libs/stats/src/uqFiniteDistribution.C \
	$(top_srcdir)/libs/stats/src/uqMetropolisHastingsSGOptions.C \
	$(top_srcdir)/libs/stats/src/uqMLSampling.C \
	$(top_srcdir)/libs/stats/src/uqMLSamplingOptions.C \
	$(top_srcdir)/libs/stats/src/uqMLSamplingLevelOptions.C \
	$(top_srcdir)/libs/stats/src/uqMonteCarloSGOptions.C \
	$(top_srcdir)/libs/stats/src/uqStatisticalInverseProblemOptions.C \
	$(top_srcdir)/libs/stats/src/uqStatisticalForwardProblemOptions.C

# ???????????????????????
libqueso_includedir = $(prefix)/include

# (5/5) Sources from interface 
libqueso_la_SOURCES += $(top_srcdir)/libs/interface/basic_classes.cpp \
	$(top_srcdir)/libs/interface/basic_interface.cpp

#libqueso_a_LIBADD   = $(HPCT_LIBSTATIC)

# (1/5) Headers to install from core/inc
libqueso_include_HEADERS = $(top_srcdir)/libs/core/inc/uqDefines.h \
	$(top_srcdir)/libs/core/inc/uqEnvironment.h \
	$(top_srcdir)/libs/core/inc/uqEnvironmentOptions.h \
	$(top_srcdir)/libs/core/inc/uqGslMatrix.h \
	$(top_srcdir)/libs/core/inc/uqGslVector.h \
	$(top_srcdir)/libs/core/inc/uqMatrix.h \
	$(top_srcdir)/libs/core/inc/uqMPI.h \
	$(top_srcdir)/libs/core/inc/uqTrilinosMatrix.h \
	$(top_srcdir)/libs/core/inc/uqTrilinosVector.h \
	$(top_srcdir)/libs/core/inc/uqVector.h

# (2/5) Headers to install from misc/inc
libqueso_include_HEADERS += $(top_srcdir)/libs/misc/inc/uq1D1DFunction.h \
	$(top_srcdir)/libs/misc/inc/uq1DQuadrature.h \
	$(top_srcdir)/libs/misc/inc/uq2dArrayOfStuff.h \
	$(top_srcdir)/libs/misc/inc/uqArrayOfOneDGrids.h \
	$(top_srcdir)/libs/misc/inc/uqArrayOfOneDTables.h \
	$(top_srcdir)/libs/misc/inc/uqAsciiTable.h \
	$(top_srcdir)/libs/misc/inc/uqCovCond.h \
	$(top_srcdir)/libs/misc/inc/uqFft.h \
	$(top_srcdir)/libs/misc/inc/uqMiscellaneous.h \
	$(top_srcdir)/libs/misc/inc/uqOneDGrid.h

# (3/5) Headers to install from basic/inc
libqueso_include_HEADERS += $(top_srcdir)/libs/basic/inc/uqArrayOfSequences.h \
	$(top_srcdir)/libs/basic/inc/uqInstantiateIntersection.h \
	$(top_srcdir)/libs/basic/inc/uqScalarFunction.h \
	$(top_srcdir)/libs/basic/inc/uqScalarFunctionSynchronizer.h \
	$(top_srcdir)/libs/basic/inc/uqScalarSequence.h \
	$(top_srcdir)/libs/basic/inc/uqSequenceOfVectors.h \
	$(top_srcdir)/libs/basic/inc/uqSequenceStatisticalOptions.h \
	$(top_srcdir)/libs/basic/inc/uqVectorFunction.h \
	$(top_srcdir)/libs/basic/inc/uqVectorFunctionSynchronizer.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSequence.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSet.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSpace.h \
	$(top_srcdir)/libs/basic/inc/uqVectorSubset.h

# (4/5) Headers to install from stats/inc
libqueso_include_HEADERS += $(top_srcdir)/libs/stats/inc/uqFiniteDistribution.h \
	$(top_srcdir)/libs/stats/inc/uqJointPdf.h \
	$(top_srcdir)/libs/stats/inc/uqMarkovChainPositionData.h \
	$(top_srcdir)/libs/stats/inc/uqMetropolisHastingsSG1.h \
	$(top_srcdir)/libs/stats/inc/uqMetropolisHastingsSG2.h \
	$(top_srcdir)/libs/stats/inc/uqMetropolisHastingsSGOptions.h \
	$(top_srcdir)/libs/stats/inc/uqMLSampling.h \
	$(top_srcdir)/libs/stats/inc/uqMLSamplingOptions.h \
	$(top_srcdir)/libs/stats/inc/uqMLSamplingLevelOptions.h \
	$(top_srcdir)/libs/stats/inc/uqModelValidation.h \
	$(top_srcdir)/libs/stats/inc/uqMonteCarloSG.h \
	$(top_srcdir)/libs/stats/inc/uqMonteCarloSGOptions.h \
	$(top_srcdir)/libs/stats/inc/uqScalarCdf.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalForwardProblem.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalForwardProblemOptions.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalInverseProblem.h \
	$(top_srcdir)/libs/stats/inc/uqStatisticalInverseProblemOptions.h \
	$(top_srcdir)/libs/stats/inc/uqTKGroup.h \
	$(top_srcdir)/libs/stats/inc/uqValidationCycle.h \
	$(top_srcdir)/libs/stats/inc/uqVectorCdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorMdf.h \
	$(top_srcdir)/libs/stats/inc/uqVectorRealizer.h \
	$(top_srcdir)/libs/stats/inc/uqVectorRV.h

# (5/5) Headers to install from interface
libqueso_include_HEADERS += $(top_srcdir)/libs/interface/queso.h \
			    $(top_srcdir)/libs/interface/basic_classes.h \
                            $(top_srcdir)/libs/interface/basic_int.h


