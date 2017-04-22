AM_CPPFLAGS = $(QUESO_CPPFLAGS)
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GSL_CFLAGS)

if HAVE_ANN
AM_CPPFLAGS += $(ANN_CFLAGS)
endif

if GRVY_ENABLED
  AM_CPPFLAGS += $(GRVY_CFLAGS)
endif

if GLPK_ENABLED
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED
  AM_CPPFLAGS += $(HDF5_CFLAGS)
endif

if LIBMESH_SLEPC_ENABLED
  AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)# $(LIBMESH_INCLUDE)
endif

BUILT_SOURCES =

# Dynamic library build: we include external library linkage here to
# make it easier for the user.

lib_LTLIBRARIES      = libqueso.la
libqueso_includedir  = $(prefix)/include/queso
libqueso_la_LDFLAGS  = $(all_libraries) -release $(GENERIC_RELEASE)
libqueso_la_LDFLAGS += $(GSL_LIBS)
libqueso_la_LDFLAGS += $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
libqueso_la_LDFLAGS += $(HDF5_LIBS)
libqueso_la_LIBADD   =  # Empty.  Append below.

if HAVE_ANN
libqueso_la_LIBADD += contrib/ANN/libANN.la
endif

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

if LIBMESH_SLEPC_ENABLED
  libqueso_la_LDFLAGS += $(LIBMESH_LIBS)
endif

# Sources from core/src

libqueso_la_SOURCES =
libqueso_la_SOURCES += core/src/Defines.C
libqueso_la_SOURCES += core/src/Environment.C
libqueso_la_SOURCES += core/src/EnvironmentOptions.C
libqueso_la_SOURCES += core/src/RngBase.C
libqueso_la_SOURCES += core/src/RngGsl.C
libqueso_la_SOURCES += core/src/RngBoost.C
libqueso_la_SOURCES += core/src/RngCXX11.C
libqueso_la_SOURCES += core/src/BasicPdfsBase.C
libqueso_la_SOURCES += core/src/BasicPdfsGsl.C
libqueso_la_SOURCES += core/src/BasicPdfsBoost.C
libqueso_la_SOURCES += core/src/BasicPdfsCXX11.C
libqueso_la_SOURCES += core/src/Vector.C
libqueso_la_SOURCES += core/src/Matrix.C
libqueso_la_SOURCES += core/src/TeuchosVector.C
libqueso_la_SOURCES += core/src/TeuchosMatrix.C
libqueso_la_SOURCES += core/src/MpiComm.C
libqueso_la_SOURCES += core/src/Map.C
libqueso_la_SOURCES += core/src/DistArray.C
libqueso_la_SOURCES += core/src/Optimizer.C
libqueso_la_SOURCES += core/src/GslOptimizer.C
libqueso_la_SOURCES += core/src/OptimizerMonitor.C
libqueso_la_SOURCES += core/src/OptimizerOptions.C
libqueso_la_SOURCES += core/src/BaseInputOptionsParser.C
libqueso_la_SOURCES += core/src/BoostInputOptionsParser.C
libqueso_la_SOURCES += core/src/TransitionKernelFactory.C
libqueso_la_SOURCES += core/src/AlgorithmFactory.C
libqueso_la_SOURCES += core/src/TKFactoryInitializer.C
libqueso_la_SOURCES += core/src/AlgorithmFactoryInitializer.C


# Sources that need libmesh
libqueso_la_SOURCES += core/src/FunctionBase.C
libqueso_la_SOURCES += core/src/LibMeshFunction.C
libqueso_la_SOURCES += core/src/OperatorBase.C
libqueso_la_SOURCES += core/src/LibMeshOperatorBase.C
libqueso_la_SOURCES += core/src/LibMeshNegativeLaplacianOperator.C
libqueso_la_SOURCES += core/src/InfiniteDimensionalGaussian.C
libqueso_la_SOURCES += core/src/InfiniteDimensionalMeasureBase.C
libqueso_la_SOURCES += core/src/InfiniteDimensionalMCMCSampler.C
libqueso_la_SOURCES += core/src/InfiniteDimensionalMCMCSamplerOptions.C
libqueso_la_SOURCES += core/src/InfiniteDimensionalLikelihoodBase.C
libqueso_la_SOURCES += core/src/FunctionOperatorBuilder.C
libqueso_la_SOURCES += core/src/GslBlockMatrix.C


# Sources from core/src with gsl conditional

if UQBT_GSL
libqueso_la_SOURCES += core/src/GslVector.C
libqueso_la_SOURCES += core/src/GslMatrix.C
endif

# Sources from misc/src

libqueso_la_SOURCES += misc/src/1D1DFunction.C
libqueso_la_SOURCES += misc/src/1DQuadrature.C
libqueso_la_SOURCES += misc/src/ComplexFft.C
libqueso_la_SOURCES += misc/src/Miscellaneous.C
libqueso_la_SOURCES += misc/src/RealFft.C
libqueso_la_SOURCES += misc/src/2dArrayOfStuff.C
libqueso_la_SOURCES += misc/src/OneDGrid.C
libqueso_la_SOURCES += misc/src/UniformOneDGrid.C
libqueso_la_SOURCES += misc/src/StdOneDGrid.C
libqueso_la_SOURCES += misc/src/Fft.C
libqueso_la_SOURCES += misc/src/CovCond.C
libqueso_la_SOURCES += misc/src/AsciiTable.C
libqueso_la_SOURCES += misc/src/ArrayOfOneDTables.C
libqueso_la_SOURCES += misc/src/ArrayOfOneDGrids.C
libqueso_la_SOURCES += misc/src/StreamUtilities.C
libqueso_la_SOURCES += misc/src/MonteCarloQuadrature.C
libqueso_la_SOURCES += misc/src/MultiDimensionalIndexing.C
libqueso_la_SOURCES += misc/src/TensorProductQuadrature.C

# Sources from misc/src withn gsl conditional

if UQBT_GSL
  libqueso_la_SOURCES += misc/src/GslAsciiTable.C
endif

# Sources from basic/src

libqueso_la_SOURCES += basic/src/TeuchosVectorSpace.C
libqueso_la_SOURCES += basic/src/SequenceStatisticalOptions.C
libqueso_la_SOURCES += basic/src/ArrayOfSequences.C
libqueso_la_SOURCES += basic/src/VectorSpace.C
libqueso_la_SOURCES += basic/src/VectorSet.C
libqueso_la_SOURCES += basic/src/BoxSubset.C
libqueso_la_SOURCES += basic/src/ConcatenationSubset.C
libqueso_la_SOURCES += basic/src/DiscreteSubset.C
libqueso_la_SOURCES += basic/src/IntersectionSubset.C
libqueso_la_SOURCES += basic/src/VectorSubset.C
libqueso_la_SOURCES += basic/src/ScalarFunction.C
libqueso_la_SOURCES += basic/src/GenericScalarFunction.C
libqueso_la_SOURCES += basic/src/ConstantScalarFunction.C
libqueso_la_SOURCES += basic/src/ScalarFunctionSynchronizer.C
libqueso_la_SOURCES += basic/src/InstantiateIntersection.C
libqueso_la_SOURCES += basic/src/SequenceOfVectors.C
libqueso_la_SOURCES += basic/src/VectorFunction.C
libqueso_la_SOURCES += basic/src/GenericVectorFunction.C
libqueso_la_SOURCES += basic/src/ConstantVectorFunction.C
libqueso_la_SOURCES += basic/src/ScalarSequence.C
libqueso_la_SOURCES += basic/src/VectorFunctionSynchronizer.C
libqueso_la_SOURCES += basic/src/VectorSequence.C


# Sources from basic/src with gsl conditional

if UQBT_GSL
  libqueso_la_SOURCES += basic/src/GslVectorSpace.C
endif

# Sources from stats/src

libqueso_la_SOURCES += stats/src/FiniteDistribution.C
libqueso_la_SOURCES += stats/src/MetropolisHastingsSG.C
libqueso_la_SOURCES += stats/src/MetropolisHastingsSGOptions.C
libqueso_la_SOURCES += stats/src/MLSampling.C
libqueso_la_SOURCES += stats/src/MLSamplingOptions.C
libqueso_la_SOURCES += stats/src/MLSamplingLevelOptions.C
libqueso_la_SOURCES += stats/src/MonteCarloSG.C
libqueso_la_SOURCES += stats/src/MonteCarloSGOptions.C
libqueso_la_SOURCES += stats/src/StatisticalInverseProblemOptions.C
libqueso_la_SOURCES += stats/src/StatisticalForwardProblem.C
libqueso_la_SOURCES += stats/src/StatisticalInverseProblem.C
libqueso_la_SOURCES += stats/src/StatisticalForwardProblemOptions.C
libqueso_la_SOURCES += stats/src/JointPdf.C
libqueso_la_SOURCES += stats/src/BayesianJointPdf.C
libqueso_la_SOURCES += stats/src/BetaJointPdf.C
libqueso_la_SOURCES += stats/src/ConcatenatedJointPdf.C
libqueso_la_SOURCES += stats/src/GammaJointPdf.C
libqueso_la_SOURCES += stats/src/GaussianJointPdf.C
libqueso_la_SOURCES += stats/src/InvLogitGaussianJointPdf.C
libqueso_la_SOURCES += stats/src/GenericJointPdf.C
libqueso_la_SOURCES += stats/src/InverseGammaJointPdf.C
libqueso_la_SOURCES += stats/src/JeffreysJointPdf.C
libqueso_la_SOURCES += stats/src/LogNormalJointPdf.C
libqueso_la_SOURCES += stats/src/PoweredJointPdf.C
libqueso_la_SOURCES += stats/src/UniformJointPdf.C
libqueso_la_SOURCES += stats/src/WignerJointPdf.C
libqueso_la_SOURCES += stats/src/MarkovChainPositionData.C
libqueso_la_SOURCES += stats/src/MatrixCovarianceFunction.C
libqueso_la_SOURCES += stats/src/ExponentialMatrixCovarianceFunction.C
libqueso_la_SOURCES += stats/src/GenericMatrixCovarianceFunction.C
libqueso_la_SOURCES += stats/src/ModelValidation.C
libqueso_la_SOURCES += stats/src/ScalarCdf.C
libqueso_la_SOURCES += stats/src/SampledScalarCdf.C
libqueso_la_SOURCES += stats/src/StdScalarCdf.C
libqueso_la_SOURCES += stats/src/ScalarCovarianceFunction.C
libqueso_la_SOURCES += stats/src/ExponentialScalarCovarianceFunction.C
libqueso_la_SOURCES += stats/src/GenericScalarCovarianceFunction.C
libqueso_la_SOURCES += stats/src/ScalarGaussianRandomField.C
libqueso_la_SOURCES += stats/src/ValidationCycle.C
libqueso_la_SOURCES += stats/src/VectorGaussianRandomField.C
libqueso_la_SOURCES += stats/src/TKGroup.C
libqueso_la_SOURCES += stats/src/ScaledCovMatrixTKGroup.C
libqueso_la_SOURCES += stats/src/TransformedScaledCovMatrixTKGroup.C
libqueso_la_SOURCES += stats/src/HessianCovMatricesTKGroup.C
libqueso_la_SOURCES += stats/src/VectorCdf.C
libqueso_la_SOURCES += stats/src/GenericVectorCdf.C
libqueso_la_SOURCES += stats/src/SampledVectorCdf.C
libqueso_la_SOURCES += stats/src/GaussianVectorCdf.C
libqueso_la_SOURCES += stats/src/VectorMdf.C
libqueso_la_SOURCES += stats/src/GenericVectorMdf.C
libqueso_la_SOURCES += stats/src/GaussianVectorMdf.C
libqueso_la_SOURCES += stats/src/SampledVectorMdf.C
libqueso_la_SOURCES += stats/src/BetaVectorRealizer.C
libqueso_la_SOURCES += stats/src/ConcatenatedVectorRealizer.C
libqueso_la_SOURCES += stats/src/GammaVectorRealizer.C
libqueso_la_SOURCES += stats/src/GaussianVectorRealizer.C
libqueso_la_SOURCES += stats/src/InvLogitGaussianVectorRealizer.C
libqueso_la_SOURCES += stats/src/GenericVectorRealizer.C
libqueso_la_SOURCES += stats/src/InverseGammaVectorRealizer.C
libqueso_la_SOURCES += stats/src/JeffreysVectorRealizer.C
libqueso_la_SOURCES += stats/src/LogNormalVectorRealizer.C
libqueso_la_SOURCES += stats/src/SequentialVectorRealizer.C
libqueso_la_SOURCES += stats/src/UniformVectorRealizer.C
libqueso_la_SOURCES += stats/src/VectorRealizer.C
libqueso_la_SOURCES += stats/src/WignerVectorRealizer.C
libqueso_la_SOURCES += stats/src/VectorRV.C
libqueso_la_SOURCES += stats/src/BetaVectorRV.C
libqueso_la_SOURCES += stats/src/ConcatenatedVectorRV.C
libqueso_la_SOURCES += stats/src/GammaVectorRV.C
libqueso_la_SOURCES += stats/src/GaussianVectorRV.C
libqueso_la_SOURCES += stats/src/InvLogitGaussianVectorRV.C
libqueso_la_SOURCES += stats/src/GenericVectorRV.C
libqueso_la_SOURCES += stats/src/InverseGammaVectorRV.C
libqueso_la_SOURCES += stats/src/JeffreysVectorRV.C
libqueso_la_SOURCES += stats/src/LogNormalVectorRV.C
libqueso_la_SOURCES += stats/src/UniformVectorRV.C
libqueso_la_SOURCES += stats/src/WignerVectorRV.C
libqueso_la_SOURCES += stats/src/LikelihoodBase.C
libqueso_la_SOURCES += stats/src/GaussianLikelihoodScalarCovariance.C
libqueso_la_SOURCES += stats/src/GaussianLikelihoodDiagonalCovariance.C
libqueso_la_SOURCES += stats/src/GaussianLikelihoodFullCovariance.C
libqueso_la_SOURCES += stats/src/GaussianLikelihoodFullCovarianceRandomCoefficient.C
libqueso_la_SOURCES += stats/src/GaussianLikelihoodBlockDiagonalCovariance.C
libqueso_la_SOURCES += stats/src/GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients.C
libqueso_la_SOURCES += stats/src/Algorithm.C
libqueso_la_SOURCES += stats/src/MetropolisAdjustedLangevinTK.C

# Sources from surrogates/src
libqueso_la_SOURCES += surrogates/src/InterpolationSurrogateData.C
libqueso_la_SOURCES += surrogates/src/InterpolationSurrogateDataSet.C
libqueso_la_SOURCES += surrogates/src/InterpolationSurrogateBase.C
libqueso_la_SOURCES += surrogates/src/LinearLagrangeInterpolationSurrogate.C
libqueso_la_SOURCES += surrogates/src/InterpolationSurrogateBuilder.C
libqueso_la_SOURCES += surrogates/src/InterpolationSurrogateIOBase.C
libqueso_la_SOURCES += surrogates/src/InterpolationSurrogateIOASCII.C

# Sources from gp/src
libqueso_la_SOURCES += gp/src/GPMSA.C
libqueso_la_SOURCES += gp/src/GPMSAOptions.C

# Headers to install from core/inc

libqueso_include_HEADERS =
libqueso_include_HEADERS += core/inc/queso.h
libqueso_include_HEADERS += core/inc/Defines.h
libqueso_include_HEADERS += core/inc/Environment.h
libqueso_include_HEADERS += core/inc/EnvironmentOptions.h
libqueso_include_HEADERS += core/inc/RngBase.h
libqueso_include_HEADERS += core/inc/RngGsl.h
libqueso_include_HEADERS += core/inc/RngBoost.h
libqueso_include_HEADERS += core/inc/RngCXX11.h
libqueso_include_HEADERS += core/inc/BasicPdfsBase.h
libqueso_include_HEADERS += core/inc/BasicPdfsGsl.h
libqueso_include_HEADERS += core/inc/BasicPdfsBoost.h
libqueso_include_HEADERS += core/inc/BasicPdfsCXX11.h
libqueso_include_HEADERS += core/inc/GslMatrix.h
libqueso_include_HEADERS += core/inc/GslVector.h
libqueso_include_HEADERS += core/inc/TeuchosMatrix.h
libqueso_include_HEADERS += core/inc/TeuchosVector.h
libqueso_include_HEADERS += core/inc/Matrix.h
libqueso_include_HEADERS += core/inc/MpiComm.h
libqueso_include_HEADERS += core/inc/Map.h
libqueso_include_HEADERS += core/inc/DistArray.h
libqueso_include_HEADERS += core/inc/asserts.h
libqueso_include_HEADERS += core/inc/exceptions.h
libqueso_include_HEADERS += core/inc/Vector.h
libqueso_include_HEADERS += core/inc/Optimizer.h
libqueso_include_HEADERS += core/inc/GslOptimizer.h
libqueso_include_HEADERS += core/inc/OptimizerMonitor.h
libqueso_include_HEADERS += core/inc/OptimizerOptions.h
libqueso_include_HEADERS += core/inc/BaseInputOptionsParser.h
libqueso_include_HEADERS += core/inc/BoostInputOptionsParser.h
libqueso_include_HEADERS += core/inc/Factory.h
libqueso_include_HEADERS += core/inc/TransitionKernelFactory.h
libqueso_include_HEADERS += core/inc/TKFactoryStochasticNewton.h
libqueso_include_HEADERS += core/inc/TKFactoryRandomWalk.h
libqueso_include_HEADERS += core/inc/TKFactoryLogitRandomWalk.h
libqueso_include_HEADERS += core/inc/AlgorithmFactory.h
libqueso_include_HEADERS += core/inc/TKFactoryMALA.h
libqueso_include_HEADERS += core/inc/TKFactoryInitializer.h
libqueso_include_HEADERS += core/inc/AlgorithmFactoryInitializer.h

# Headers that need libmesh
libqueso_include_HEADERS += core/inc/FunctionBase.h
libqueso_include_HEADERS += core/inc/LibMeshFunction.h
libqueso_include_HEADERS += core/inc/OperatorBase.h
libqueso_include_HEADERS += core/inc/LibMeshOperatorBase.h
libqueso_include_HEADERS += core/inc/LibMeshNegativeLaplacianOperator.h
libqueso_include_HEADERS += core/inc/InfiniteDimensionalMeasureBase.h
libqueso_include_HEADERS += core/inc/InfiniteDimensionalGaussian.h
libqueso_include_HEADERS += core/inc/InfiniteDimensionalMCMCSampler.h
libqueso_include_HEADERS += core/inc/InfiniteDimensionalMCMCSamplerOptions.h
libqueso_include_HEADERS += core/inc/InfiniteDimensionalLikelihoodBase.h
libqueso_include_HEADERS += core/inc/FunctionOperatorBuilder.h
libqueso_include_HEADERS += core/inc/GslBlockMatrix.h
libqueso_include_HEADERS += core/inc/ScopedPtr.h
libqueso_include_HEADERS += core/inc/SharedPtr.h

# Headers to install from misc/inc

libqueso_include_HEADERS += misc/inc/BaseQuadrature.h
libqueso_include_HEADERS += misc/inc/1D1DFunction.h
libqueso_include_HEADERS += misc/inc/1DQuadrature.h
libqueso_include_HEADERS += misc/inc/2dArrayOfStuff.h
libqueso_include_HEADERS += misc/inc/ArrayOfOneDGrids.h
libqueso_include_HEADERS += misc/inc/ArrayOfOneDTables.h
libqueso_include_HEADERS += misc/inc/AsciiTable.h
libqueso_include_HEADERS += misc/inc/CovCond.h
libqueso_include_HEADERS += misc/inc/Fft.h
libqueso_include_HEADERS += misc/inc/Miscellaneous.h
libqueso_include_HEADERS += misc/inc/OneDGrid.h
libqueso_include_HEADERS += misc/inc/UniformOneDGrid.h
libqueso_include_HEADERS += misc/inc/StdOneDGrid.h
libqueso_include_HEADERS += misc/inc/StreamUtilities.h
libqueso_include_HEADERS += misc/inc/math_macros.h
libqueso_include_HEADERS += misc/inc/MultiDQuadratureBase.h
libqueso_include_HEADERS += misc/inc/MonteCarloQuadrature.h
libqueso_include_HEADERS += misc/inc/MultiDimensionalIndexing.h
libqueso_include_HEADERS += misc/inc/TensorProductQuadrature.h

# Headers to install from basic/inc

libqueso_include_HEADERS += basic/inc/ArrayOfSequences.h
libqueso_include_HEADERS += basic/inc/InstantiateIntersection.h
libqueso_include_HEADERS += basic/inc/ScalarFunction.h
libqueso_include_HEADERS += basic/inc/GenericScalarFunction.h
libqueso_include_HEADERS += basic/inc/ConstantScalarFunction.h
libqueso_include_HEADERS += basic/inc/ScalarFunctionSynchronizer.h
libqueso_include_HEADERS += basic/inc/ScalarSequence.h
libqueso_include_HEADERS += basic/inc/SequenceOfVectors.h
libqueso_include_HEADERS += basic/inc/SequenceStatisticalOptions.h
libqueso_include_HEADERS += basic/inc/VectorFunction.h
libqueso_include_HEADERS += basic/inc/VectorFunctionSynchronizer.h
libqueso_include_HEADERS += basic/inc/VectorSequence.h
libqueso_include_HEADERS += basic/inc/VectorSet.h
libqueso_include_HEADERS += basic/inc/VectorSpace.h
libqueso_include_HEADERS += basic/inc/VectorSubset.h
libqueso_include_HEADERS += basic/inc/BoxSubset.h
libqueso_include_HEADERS += basic/inc/ConcatenationSubset.h
libqueso_include_HEADERS += basic/inc/DiscreteSubset.h
libqueso_include_HEADERS += basic/inc/IntersectionSubset.h
libqueso_include_HEADERS += basic/inc/GenericVectorFunction.h
libqueso_include_HEADERS += basic/inc/ConstantVectorFunction.h

# Headers to install from stats/inc

libqueso_include_HEADERS += stats/inc/FiniteDistribution.h
libqueso_include_HEADERS += stats/inc/JointPdf.h
libqueso_include_HEADERS += stats/inc/BayesianJointPdf.h
libqueso_include_HEADERS += stats/inc/BetaJointPdf.h
libqueso_include_HEADERS += stats/inc/ConcatenatedJointPdf.h
libqueso_include_HEADERS += stats/inc/GammaJointPdf.h
libqueso_include_HEADERS += stats/inc/GaussianJointPdf.h
libqueso_include_HEADERS += stats/inc/InvLogitGaussianJointPdf.h
libqueso_include_HEADERS += stats/inc/GenericJointPdf.h
libqueso_include_HEADERS += stats/inc/InverseGammaJointPdf.h
libqueso_include_HEADERS += stats/inc/JeffreysJointPdf.h
libqueso_include_HEADERS += stats/inc/LogNormalJointPdf.h
libqueso_include_HEADERS += stats/inc/PoweredJointPdf.h
libqueso_include_HEADERS += stats/inc/UniformJointPdf.h
libqueso_include_HEADERS += stats/inc/WignerJointPdf.h
libqueso_include_HEADERS += stats/inc/MarkovChainPositionData.h
libqueso_include_HEADERS += stats/inc/MetropolisHastingsSG.h
libqueso_include_HEADERS += stats/inc/MetropolisHastingsSGOptions.h
libqueso_include_HEADERS += stats/inc/MLSampling.h
libqueso_include_HEADERS += stats/inc/MLSamplingOptions.h
libqueso_include_HEADERS += stats/inc/MLSamplingLevelOptions.h
libqueso_include_HEADERS += stats/inc/ModelValidation.h
libqueso_include_HEADERS += stats/inc/MonteCarloSG.h
libqueso_include_HEADERS += stats/inc/MonteCarloSGOptions.h
libqueso_include_HEADERS += stats/inc/ScalarCdf.h
libqueso_include_HEADERS += stats/inc/SampledScalarCdf.h
libqueso_include_HEADERS += stats/inc/StdScalarCdf.h
libqueso_include_HEADERS += stats/inc/StatisticalForwardProblem.h
libqueso_include_HEADERS += stats/inc/StatisticalForwardProblemOptions.h
libqueso_include_HEADERS += stats/inc/StatisticalInverseProblem.h
libqueso_include_HEADERS += stats/inc/StatisticalInverseProblemOptions.h
libqueso_include_HEADERS += stats/inc/TKGroup.h
libqueso_include_HEADERS += stats/inc/ScaledCovMatrixTKGroup.h
libqueso_include_HEADERS += stats/inc/TransformedScaledCovMatrixTKGroup.h
libqueso_include_HEADERS += stats/inc/HessianCovMatricesTKGroup.h
libqueso_include_HEADERS += stats/inc/ValidationCycle.h
libqueso_include_HEADERS += stats/inc/VectorCdf.h
libqueso_include_HEADERS += stats/inc/GenericVectorCdf.h
libqueso_include_HEADERS += stats/inc/SampledVectorCdf.h
libqueso_include_HEADERS += stats/inc/GaussianVectorCdf.h
libqueso_include_HEADERS += stats/inc/VectorMdf.h
libqueso_include_HEADERS += stats/inc/GenericVectorMdf.h
libqueso_include_HEADERS += stats/inc/GaussianVectorMdf.h
libqueso_include_HEADERS += stats/inc/SampledVectorMdf.h
libqueso_include_HEADERS += stats/inc/BetaVectorRealizer.h
libqueso_include_HEADERS += stats/inc/ConcatenatedVectorRealizer.h
libqueso_include_HEADERS += stats/inc/GammaVectorRealizer.h
libqueso_include_HEADERS += stats/inc/GaussianVectorRealizer.h
libqueso_include_HEADERS += stats/inc/InvLogitGaussianVectorRealizer.h
libqueso_include_HEADERS += stats/inc/GenericVectorRealizer.h
libqueso_include_HEADERS += stats/inc/InverseGammaVectorRealizer.h
libqueso_include_HEADERS += stats/inc/JeffreysVectorRealizer.h
libqueso_include_HEADERS += stats/inc/LogNormalVectorRealizer.h
libqueso_include_HEADERS += stats/inc/SequentialVectorRealizer.h
libqueso_include_HEADERS += stats/inc/UniformVectorRealizer.h
libqueso_include_HEADERS += stats/inc/VectorRealizer.h
libqueso_include_HEADERS += stats/inc/WignerVectorRealizer.h
libqueso_include_HEADERS += stats/inc/VectorRV.h
libqueso_include_HEADERS += stats/inc/BetaVectorRV.h
libqueso_include_HEADERS += stats/inc/ConcatenatedVectorRV.h
libqueso_include_HEADERS += stats/inc/GammaVectorRV.h
libqueso_include_HEADERS += stats/inc/GaussianVectorRV.h
libqueso_include_HEADERS += stats/inc/InvLogitGaussianVectorRV.h
libqueso_include_HEADERS += stats/inc/GenericVectorRV.h
libqueso_include_HEADERS += stats/inc/InverseGammaVectorRV.h
libqueso_include_HEADERS += stats/inc/JeffreysVectorRV.h
libqueso_include_HEADERS += stats/inc/LogNormalVectorRV.h
libqueso_include_HEADERS += stats/inc/UniformVectorRV.h
libqueso_include_HEADERS += stats/inc/WignerVectorRV.h
libqueso_include_HEADERS += stats/inc/InfoTheory.h
libqueso_include_HEADERS += stats/inc/InfoTheory_impl.h
libqueso_include_HEADERS += stats/inc/ScalarCovarianceFunction.h
libqueso_include_HEADERS += stats/inc/ExponentialScalarCovarianceFunction.h
libqueso_include_HEADERS += stats/inc/GenericScalarCovarianceFunction.h
libqueso_include_HEADERS += stats/inc/MatrixCovarianceFunction.h
libqueso_include_HEADERS += stats/inc/ExponentialMatrixCovarianceFunction.h
libqueso_include_HEADERS += stats/inc/GenericMatrixCovarianceFunction.h
libqueso_include_HEADERS += stats/inc/ScalarGaussianRandomField.h
libqueso_include_HEADERS += stats/inc/VectorGaussianRandomField.h
libqueso_include_HEADERS += stats/inc/LikelihoodBase.h
libqueso_include_HEADERS += stats/inc/GaussianLikelihoodScalarCovariance.h
libqueso_include_HEADERS += stats/inc/GaussianLikelihoodDiagonalCovariance.h
libqueso_include_HEADERS += stats/inc/GaussianLikelihoodFullCovariance.h
libqueso_include_HEADERS += stats/inc/GaussianLikelihoodFullCovarianceRandomCoefficient.h
libqueso_include_HEADERS += stats/inc/GaussianLikelihoodBlockDiagonalCovariance.h
libqueso_include_HEADERS += stats/inc/GaussianLikelihoodBlockDiagonalCovarianceRandomCoefficients.h
libqueso_include_HEADERS += stats/inc/Algorithm.h
libqueso_include_HEADERS += stats/inc/MetropolisAdjustedLangevinTK.h

# Headers to install from surrogates/inc
libqueso_include_HEADERS += surrogates/inc/SurrogateBase.h
libqueso_include_HEADERS += surrogates/inc/InterpolationSurrogateData.h
libqueso_include_HEADERS += surrogates/inc/InterpolationSurrogateDataSet.h
libqueso_include_HEADERS += surrogates/inc/InterpolationSurrogateBase.h
libqueso_include_HEADERS += surrogates/inc/LinearLagrangeInterpolationSurrogate.h
libqueso_include_HEADERS += surrogates/inc/SurrogateBuilderBase.h
libqueso_include_HEADERS += surrogates/inc/InterpolationSurrogateBuilder.h
libqueso_include_HEADERS += surrogates/inc/InterpolationSurrogateIOBase.h
libqueso_include_HEADERS += surrogates/inc/InterpolationSurrogateIOASCII.h

# Headers to install from gp/inc
libqueso_include_HEADERS += gp/inc/GPMSA.h
libqueso_include_HEADERS += gp/inc/GPMSAOptions.h

# Install the magic header
libqueso_include_HEADERS += contrib/inc/all.h
libqueso_include_HEADERS += contrib/getpot/getpot.h

# Headers to install, but not part of dist

nodist_libqueso_include_HEADERS =
nodist_libqueso_include_HEADERS += $(top_builddir)/config_queso.h

# Standalone binaries

bin_PROGRAMS          = queso_version
queso_version_SOURCES = core/src/version.C
queso_version_LDADD   = libqueso.la

if TRILINOS_ENABLED
  queso_version_LDADD     += -lteuchoscore -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosremainder -lepetra
endif

queso_version_DEPENDENCIES = libqueso.la

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif
