# Build in these directories:
# DIST_SUBDIRS = basic stats

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += -I$(top_builddir)/src/core/inc  # For queso.h
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CFLAGS)
AM_CPPFLAGS += $(ANN_CFLAGS)

###	      -I$(top_srcdir)/src/interface 

if GRVY_ENABLED
  AM_CPPFLAGS += $(GRVY_CFLAGS)
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

if LIBMESH_ENABLED
  AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)# $(LIBMESH_INCLUDE)
endif

if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

# Dynamic library build: we include external library linkage here to
# make it easier for the user.

lib_LTLIBRARIES      = libqueso.la
libqueso_includedir  = $(prefix)/include/queso
libqueso_la_LDFLAGS  = $(all_libraries) -release $(GENERIC_RELEASE)
libqueso_la_LDFLAGS += $(GSL_LIBS) 
libqueso_la_LDFLAGS += $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
libqueso_la_LDFLAGS += $(ANN_LIBS)

if GRVY_ENABLED	
  libqueso_la_LDFLAGS += $(GRVY_LIBS)
endif

if GLPK_ENABLED	
  libqueso_la_LDFLAGS += $(GLPK_LIBS)
endif

if HDF5_ENABLED	
  libqueso_la_LDFLAGS += $(HDF5_LIBS)
endif

# TODO: cleanup way epetra is handled; it looks like it is being put
# in LDFLAGS directly. Would like to have a EPETRA_LIBS or some such
# variable later.

if TRILINOS_ENABLED
  libqueso_la_LDFLAGS += -lteuchoscore -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosremainder -lepetra
endif

if LIBMESH_ENABLED
  libqueso_la_LDFLAGS += $(LIBMESH_LIBS)
endif

#libqueso_la_LDFLAGS += -Wl,-rpath,$(TRILINOS_HOME)/lib -lepetra

# Sources from core/src

libqueso_la_SOURCES = \
	$(top_srcdir)/src/core/src/Defines.C \
	$(top_srcdir)/src/core/src/Environment.C \
	$(top_srcdir)/src/core/src/EnvironmentOptions.C \
	$(top_srcdir)/src/core/src/RngBase.C \
	$(top_srcdir)/src/core/src/RngGsl.C \
	$(top_srcdir)/src/core/src/RngBoost.C \
	$(top_srcdir)/src/core/src/BasicPdfsBase.C \
	$(top_srcdir)/src/core/src/BasicPdfsGsl.C \
	$(top_srcdir)/src/core/src/BasicPdfsBoost.C \
	$(top_srcdir)/src/core/src/Vector.C \
	$(top_srcdir)/src/core/src/Matrix.C \
	$(top_srcdir)/src/core/src/TeuchosVector.C \
	$(top_srcdir)/src/core/src/TeuchosMatrix.C \
	$(top_srcdir)/src/core/src/MpiComm.C \
	$(top_srcdir)/src/core/src/Map.C \
	$(top_srcdir)/src/core/src/DistArray.C \
	$(top_srcdir)/src/core/src/uqFunctionBase.C \
	$(top_srcdir)/src/core/src/uqLibMeshFunction.C \
	$(top_srcdir)/src/core/src/uqOperatorBase.C \
	$(top_srcdir)/src/core/src/uqLibMeshOperatorBase.C \
	$(top_srcdir)/src/core/src/uqLibMeshNegativeLaplacianOperator.C


# Sources from core/src with gsl conditional 

if UQBT_GSL
libqueso_la_SOURCES += \
	$(top_srcdir)/src/core/src/GslVector.C \
	$(top_srcdir)/src/core/src/GslMatrix.C
endif

# Sources from misc/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/misc/src/1D1DFunction.C \
	$(top_srcdir)/src/misc/src/1DQuadrature.C \
	$(top_srcdir)/src/misc/src/ComplexFft.C \
	$(top_srcdir)/src/misc/src/Miscellaneous.C \
	$(top_srcdir)/src/misc/src/RealFft.C \
	$(top_srcdir)/src/misc/src/2dArrayOfStuff.C

# Sources from misc/src withn gsl conditional

if UQBT_GSL
  libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/GslAsciiTable.C
endif

# Sources from basic/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/basic/src/TeuchosVectorSpace.C \
	$(top_srcdir)/src/basic/src/SequenceStatisticalOptions.C \
	$(top_srcdir)/src/basic/src/ArrayOfSequences.C \
	$(top_srcdir)/src/basic/src/VectorSpace.C \
	$(top_srcdir)/src/basic/src/VectorSet.C \
	$(top_srcdir)/src/basic/src/BoxSubset.C \
	$(top_srcdir)/src/basic/src/ConcatenationSubset.C \
	$(top_srcdir)/src/basic/src/DiscreteSubset.C \
	$(top_srcdir)/src/basic/src/IntersectionSubset.C \
	$(top_srcdir)/src/basic/src/VectorSubset.C


# Sources from basic/src with gsl conditional

if UQBT_GSL
  libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/GslVectorSpace.C
endif

# Sources from stats/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/stats/src/FiniteDistribution.C \
	$(top_srcdir)/src/stats/src/MetropolisHastingsSG.C \
	$(top_srcdir)/src/stats/src/MetropolisHastingsSGOptions.C \
	$(top_srcdir)/src/stats/src/MLSampling.C \
	$(top_srcdir)/src/stats/src/MLSamplingOptions.C \
	$(top_srcdir)/src/stats/src/MLSamplingLevelOptions.C \
	$(top_srcdir)/src/stats/src/MonteCarloSGOptions.C \
	$(top_srcdir)/src/stats/src/StatisticalInverseProblemOptions.C \
	$(top_srcdir)/src/stats/src/StatisticalForwardProblemOptions.C \
	$(top_srcdir)/src/stats/src/InfoTheory.C

# Sources from gp/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/gp/src/ExperimentModelOptions.C \
	$(top_srcdir)/src/gp/src/GpmsaComputerModelOptions.C \
	$(top_srcdir)/src/gp/src/SimulationModelOptions.C

# Headers to install from core/inc

libqueso_include_HEADERS = \
	$(top_builddir)/src/core/inc/queso.h \
	$(top_srcdir)/src/core/inc/Defines.h \
	$(top_srcdir)/src/core/inc/Environment.h \
	$(top_srcdir)/src/core/inc/EnvironmentOptions.h \
	$(top_srcdir)/src/core/inc/RngBase.h \
	$(top_srcdir)/src/core/inc/RngGsl.h \
	$(top_srcdir)/src/core/inc/RngBoost.h \
	$(top_srcdir)/src/core/inc/BasicPdfsBase.h \
	$(top_srcdir)/src/core/inc/BasicPdfsGsl.h \
	$(top_srcdir)/src/core/inc/BasicPdfsBoost.h \
	$(top_srcdir)/src/core/inc/GslMatrix.h \
	$(top_srcdir)/src/core/inc/GslVector.h \
	$(top_srcdir)/src/core/inc/TeuchosMatrix.h \
	$(top_srcdir)/src/core/inc/TeuchosVector.h \
	$(top_srcdir)/src/core/inc/Matrix.h \
	$(top_srcdir)/src/core/inc/MpiComm.h \
	$(top_srcdir)/src/core/inc/Map.h \
	$(top_srcdir)/src/core/inc/DistArray.h \
	$(top_srcdir)/src/core/inc/asserts.h \
	$(top_srcdir)/src/core/inc/exceptions.h \
	$(top_srcdir)/src/core/inc/Vector.h \
	$(top_srcdir)/src/core/inc/uqFunctionBase.h \
	$(top_srcdir)/src/core/inc/uqLibMeshFunction.h \
	$(top_srcdir)/src/core/inc/uqOperatorBase.h \
	$(top_srcdir)/src/core/inc/uqLibMeshOperatorBase.h \
	$(top_srcdir)/src/core/inc/uqLibMeshNegativeLaplacianOperator.h

# Headers to install from misc/inc

libqueso_include_HEADERS += \
	$(top_srcdir)/src/misc/inc/1D1DFunction.h \
	$(top_srcdir)/src/misc/inc/1DQuadrature.h \
	$(top_srcdir)/src/misc/inc/2dArrayOfStuff.h \
	$(top_srcdir)/src/misc/inc/ArrayOfOneDGrids.h \
	$(top_srcdir)/src/misc/inc/ArrayOfOneDTables.h \
	$(top_srcdir)/src/misc/inc/AsciiTable.h \
	$(top_srcdir)/src/misc/inc/CovCond.h \
	$(top_srcdir)/src/misc/inc/Fft.h \
	$(top_srcdir)/src/misc/inc/Miscellaneous.h \
	$(top_srcdir)/src/misc/inc/OneDGrid.h

# Headers to install from basic/inc

libqueso_include_HEADERS += \
	$(top_srcdir)/src/basic/inc/ArrayOfSequences.h \
	$(top_srcdir)/src/basic/inc/InstantiateIntersection.h \
	$(top_srcdir)/src/basic/inc/ScalarFunction.h \
	$(top_srcdir)/src/basic/inc/ScalarFunctionSynchronizer.h \
	$(top_srcdir)/src/basic/inc/ScalarSequence.h \
	$(top_srcdir)/src/basic/inc/SequenceOfVectors.h \
	$(top_srcdir)/src/basic/inc/SequenceStatisticalOptions.h \
	$(top_srcdir)/src/basic/inc/VectorFunction.h \
	$(top_srcdir)/src/basic/inc/VectorFunctionSynchronizer.h \
	$(top_srcdir)/src/basic/inc/VectorSequence.h \
	$(top_srcdir)/src/basic/inc/VectorSet.h \
	$(top_srcdir)/src/basic/inc/VectorSpace.h \
	$(top_srcdir)/src/basic/inc/VectorSubset.h \
	$(top_srcdir)/src/basic/inc/BoxSubset.h \
	$(top_srcdir)/src/basic/inc/ConcatenationSubset.h \
	$(top_srcdir)/src/basic/inc/DiscreteSubset.h \
	$(top_srcdir)/src/basic/inc/IntersectionSubset.h

# Headers to install from stats/inc 

libqueso_include_HEADERS += \
	$(top_srcdir)/src/stats/inc/FiniteDistribution.h \
	$(top_srcdir)/src/stats/inc/JointPdf.h \
	$(top_srcdir)/src/stats/inc/MarkovChainPositionData.h \
	$(top_srcdir)/src/stats/inc/MetropolisHastingsSG1.h \
	$(top_srcdir)/src/stats/inc/MetropolisHastingsSG2.h \
	$(top_srcdir)/src/stats/inc/MetropolisHastingsSGOptions.h \
	$(top_srcdir)/src/stats/inc/MLSampling1.h \
	$(top_srcdir)/src/stats/inc/MLSampling2.h \
	$(top_srcdir)/src/stats/inc/MLSampling3.h \
	$(top_srcdir)/src/stats/inc/MLSamplingOptions.h \
	$(top_srcdir)/src/stats/inc/MLSamplingLevelOptions.h \
	$(top_srcdir)/src/stats/inc/ModelValidation.h \
	$(top_srcdir)/src/stats/inc/MonteCarloSG.h \
	$(top_srcdir)/src/stats/inc/MonteCarloSGOptions.h \
	$(top_srcdir)/src/stats/inc/ScalarCdf.h \
	$(top_srcdir)/src/stats/inc/StatisticalForwardProblem.h \
	$(top_srcdir)/src/stats/inc/StatisticalForwardProblemOptions.h \
	$(top_srcdir)/src/stats/inc/StatisticalInverseProblem.h \
	$(top_srcdir)/src/stats/inc/StatisticalInverseProblemOptions.h \
	$(top_srcdir)/src/stats/inc/TKGroup.h \
	$(top_srcdir)/src/stats/inc/ValidationCycle.h \
	$(top_srcdir)/src/stats/inc/VectorCdf.h \
	$(top_srcdir)/src/stats/inc/VectorMdf.h \
	$(top_srcdir)/src/stats/inc/VectorRealizer.h \
	$(top_srcdir)/src/stats/inc/VectorRV.h \
	$(top_srcdir)/src/stats/inc/InfoTheory.h \
	$(top_srcdir)/src/stats/inc/ScalarCovarianceFunction.h \
	$(top_srcdir)/src/stats/inc/MatrixCovarianceFunction.h \
	$(top_srcdir)/src/stats/inc/ScalarGaussianRandomField.h \
	$(top_srcdir)/src/stats/inc/VectorGaussianRandomField.h

# Headers to install from gp/inc 

libqueso_include_HEADERS += \
	$(top_srcdir)/src/gp/inc/ExperimentModel.h \
	$(top_srcdir)/src/gp/inc/ExperimentModelOptions.h \
	$(top_srcdir)/src/gp/inc/ExperimentStorage.h \
	$(top_srcdir)/src/gp/inc/GcmExperimentInfo.h \
	$(top_srcdir)/src/gp/inc/GcmJointInfo.h \
	$(top_srcdir)/src/gp/inc/GcmJointTildeInfo.h \
	$(top_srcdir)/src/gp/inc/GcmSimulationInfo.h \
	$(top_srcdir)/src/gp/inc/GcmSimulationTildeInfo.h \
	$(top_srcdir)/src/gp/inc/GcmTotalInfo.h \
	$(top_srcdir)/src/gp/inc/GcmZInfo.h \
	$(top_srcdir)/src/gp/inc/GcmZTildeInfo.h \
	$(top_srcdir)/src/gp/inc/GpmsaComputerModel.h \
	$(top_srcdir)/src/gp/inc/GpmsaComputerModel2.h \
	$(top_srcdir)/src/gp/inc/GpmsaComputerModel3.h \
	$(top_srcdir)/src/gp/inc/GpmsaComputerModel4.h \
	$(top_srcdir)/src/gp/inc/GpmsaComputerModel5.h \
	$(top_srcdir)/src/gp/inc/GpmsaComputerModelOptions.h \
	$(top_srcdir)/src/gp/inc/SimulationModel.h \
	$(top_srcdir)/src/gp/inc/SimulationModelOptions.h \
	$(top_srcdir)/src/gp/inc/SimulationStorage.h

# Headers to install, but not part of dist

nodist_libqueso_include_HEADERS    =       \
	$(top_builddir)/config_queso.h       \
	$(top_builddir)/src/core/inc/queso.h

# Standalone binaries

bin_PROGRAMS               = queso_version
queso_version_SOURCES      = $(top_srcdir)/src/core/src/version.C
queso_version_LDADD        = -L$(top_builddir)/src -lqueso $(HDF5_LIBS) 
if TRILINOS_ENABLED
  queso_version_LDADD     += -lteuchoscore -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosremainder -lepetra
endif

queso_version_DEPENDENCIES = libqueso.la

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(libqueso_include_HEADERS) $(libqueso_la_SOURCES) $(top_srcdir)/src/core/inc/queso.h
	echo 'updated source license headers' >$@
endif
