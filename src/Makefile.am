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
libqueso_la_LDFLAGS += $(HDF5_LIBS) $(HDF5_CXXLIBS)

if GRVY_ENABLED	
  libqueso_la_LDFLAGS += $(GRVY_LIBS)
endif

if GLPK_ENABLED	
  libqueso_la_LDFLAGS += $(GLPK_LIBS)
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

libqueso_la_SOURCES =
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/Defines.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/Environment.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/EnvironmentOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/RngBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/RngGsl.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/RngBoost.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/BasicPdfsBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/BasicPdfsGsl.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/BasicPdfsBoost.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/Vector.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/Matrix.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/TeuchosVector.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/TeuchosMatrix.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/MpiComm.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/Map.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/DistArray.C


# Sources that need libmesh
if LIBMESH_ENABLED
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/FunctionBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/LibMeshFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/OperatorBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/LibMeshOperatorBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/LibMeshNegativeLaplacianOperator.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/InfiniteDimensionalGaussian.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/InfiniteDimensionalMeasureBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/InfiniteDimensionalMCMCSampler.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/InfiniteDimensionalMCMCSamplerOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/InfiniteDimensionalLikelihoodBase.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/FunctionOperatorBuilder.C
endif


# Sources from core/src with gsl conditional 

if UQBT_GSL
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/GslVector.C
libqueso_la_SOURCES += $(top_srcdir)/src/core/src/GslMatrix.C
endif

# Sources from misc/src

libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/1D1DFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/1DQuadrature.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/ComplexFft.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/Miscellaneous.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/RealFft.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/2dArrayOfStuff.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/OneDGrid.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/UniformOneDGrid.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/StdOneDGrid.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/Fft.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/CovCond.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/AsciiTable.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/ArrayOfOneDTables.C
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/ArrayOfOneDGrids.C

# Sources from misc/src withn gsl conditional

if UQBT_GSL
  libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/GslAsciiTable.C
endif

# Sources from basic/src

libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/TeuchosVectorSpace.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/SequenceStatisticalOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ArrayOfSequences.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/VectorSpace.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/VectorSet.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/BoxSubset.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ConcatenationSubset.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/DiscreteSubset.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/IntersectionSubset.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/VectorSubset.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ScalarFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/GenericScalarFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ConstantScalarFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ScalarFunctionSynchronizer.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/InstantiateIntersection.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/SequenceOfVectors.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/VectorFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/GenericVectorFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ConstantVectorFunction.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/ScalarSequence.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/VectorFunctionSynchronizer.C
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/VectorSequence.C


# Sources from basic/src with gsl conditional

if UQBT_GSL
  libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/GslVectorSpace.C
endif

# Sources from stats/src

libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/FiniteDistribution.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MetropolisHastingsSG.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MetropolisHastingsSGOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MLSampling.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MLSamplingOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MLSamplingLevelOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MonteCarloSGOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/StatisticalInverseProblemOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/StatisticalForwardProblemOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/InfoTheory.C
libqueso_la_SOURCES += $(top_srcdir)/src/stats/src/MarkovChainPositionData.C

# Sources from gp/src

libqueso_la_SOURCES += $(top_srcdir)/src/gp/src/ExperimentModelOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/gp/src/GpmsaComputerModelOptions.C
libqueso_la_SOURCES += $(top_srcdir)/src/gp/src/SimulationModelOptions.C

# Headers to install from core/inc

libqueso_include_HEADERS =
libqueso_include_HEADERS += $(top_builddir)/src/core/inc/queso.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/Defines.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/Environment.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/EnvironmentOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/RngBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/RngGsl.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/RngBoost.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/BasicPdfsBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/BasicPdfsGsl.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/BasicPdfsBoost.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/GslMatrix.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/GslVector.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/TeuchosMatrix.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/TeuchosVector.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/Matrix.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/MpiComm.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/Map.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/DistArray.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/asserts.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/exceptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/Vector.h

# Headers that need libmesh
if LIBMESH_ENABLED
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/FunctionBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/LibMeshFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/OperatorBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/LibMeshOperatorBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/LibMeshNegativeLaplacianOperator.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/InfiniteDimensionalMeasureBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/InfiniteDimensionalGaussian.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/InfiniteDimensionalMCMCSampler.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/InfiniteDimensionalMCMCSamplerOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/InfiniteDimensionalLikelihoodBase.h
libqueso_include_HEADERS += $(top_srcdir)/src/core/inc/FunctionOperatorBuilder.h
endif

# Headers to install from misc/inc

libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/1D1DFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/1DQuadrature.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/2dArrayOfStuff.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/ArrayOfOneDGrids.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/ArrayOfOneDTables.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/AsciiTable.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/CovCond.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/Fft.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/Miscellaneous.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/OneDGrid.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/UniformOneDGrid.h
libqueso_include_HEADERS += $(top_srcdir)/src/misc/inc/StdOneDGrid.h

# Headers to install from basic/inc

libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ArrayOfSequences.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/InstantiateIntersection.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ScalarFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/GenericScalarFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ConstantScalarFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ScalarFunctionSynchronizer.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ScalarSequence.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/SequenceOfVectors.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/SequenceStatisticalOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/VectorFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/VectorFunctionSynchronizer.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/VectorSequence.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/VectorSet.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/VectorSpace.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/VectorSubset.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/BoxSubset.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ConcatenationSubset.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/DiscreteSubset.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/IntersectionSubset.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/GenericVectorFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/basic/inc/ConstantVectorFunction.h

# Headers to install from stats/inc 

libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/FiniteDistribution.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/JointPdf.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MarkovChainPositionData.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MetropolisHastingsSG1.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MetropolisHastingsSG2.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MetropolisHastingsSGOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MLSampling1.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MLSampling2.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MLSampling3.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MLSamplingOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MLSamplingLevelOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/ModelValidation.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MonteCarloSG.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MonteCarloSGOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/ScalarCdf.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/StatisticalForwardProblem.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/StatisticalForwardProblemOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/StatisticalInverseProblem.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/StatisticalInverseProblemOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/TKGroup.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/ValidationCycle.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/VectorCdf.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/VectorMdf.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/VectorRealizer.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/VectorRV.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/InfoTheory.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/ScalarCovarianceFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/MatrixCovarianceFunction.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/ScalarGaussianRandomField.h
libqueso_include_HEADERS += $(top_srcdir)/src/stats/inc/VectorGaussianRandomField.h

# Headers to install from gp/inc 

libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/ExperimentModel.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/ExperimentModelOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/ExperimentStorage.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmExperimentInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmJointInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmJointTildeInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmSimulationInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmSimulationTildeInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmTotalInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmZInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GcmZTildeInfo.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GpmsaComputerModel.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GpmsaComputerModel2.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GpmsaComputerModel3.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GpmsaComputerModel4.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GpmsaComputerModel5.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/GpmsaComputerModelOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/SimulationModel.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/SimulationModelOptions.h
libqueso_include_HEADERS += $(top_srcdir)/src/gp/inc/SimulationStorage.h

# Headers to install, but not part of dist

nodist_libqueso_include_HEADERS =
nodist_libqueso_include_HEADERS += $(top_builddir)/config_queso.h
nodist_libqueso_include_HEADERS += $(top_builddir)/src/core/inc/queso.h

# Standalone binaries

bin_PROGRAMS               = queso_version
queso_version_SOURCES      = $(top_srcdir)/src/core/src/version.C
#queso_version_LDADD        = -L$(top_builddir)/src -lqueso $(HDF5_LIBS) $(HDF5_LIB)/libhdf5_cpp.la 
queso_version_LDADD = libqueso.la
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
