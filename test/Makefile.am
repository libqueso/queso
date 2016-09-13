SUBDIRS    = t01_valid_cycle t02_sip_sfp t03_sequence t04_bimodal

if GIT_CLONE
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += test_boxsubset_centroid
check_PROGRAMS += test_concatenation_centroid
check_PROGRAMS += test_beta_mean
check_PROGRAMS += test_concatenated_mean
check_PROGRAMS += test_gamma_mean
check_PROGRAMS += test_uqEnvironmentNonFatal
check_PROGRAMS += test_uqGslVector
check_PROGRAMS += test_uqGaussianVectorRVClass
check_PROGRAMS += test_VectorPdf_gsl
check_PROGRAMS += test_VectorRV_gsl
check_PROGRAMS += test_VectorRealizer_gsl
check_PROGRAMS += test_uqGslMatrix
check_PROGRAMS += test_uqTeuchosVector
check_PROGRAMS += test_uqexception
check_PROGRAMS += test_DistArrayCopy
check_PROGRAMS += test_DistArrayEquals
check_PROGRAMS += test_DistArrayMisc
check_PROGRAMS += test_GslVectorSpaceMisc
check_PROGRAMS += test_IntersectionSubsetContains
check_PROGRAMS += test_operator
check_PROGRAMS += test_inf_gaussian
check_PROGRAMS += test_inf_options
check_PROGRAMS += test_SequenceOfVectorsErase
check_PROGRAMS += test_GaussianMean1DRegression
check_PROGRAMS += test_gpmsa_cobra
check_PROGRAMS += test_gpmsa_vector
check_PROGRAMS += test_LlhdTargetOutput
check_PROGRAMS += test_jeffreys
check_PROGRAMS += test_gsloptimizer
check_PROGRAMS += test_seedwithmap
check_PROGRAMS += test_seedwithmap_fd
check_PROGRAMS += test_logitadaptedcov
check_PROGRAMS += test_scalarCovariance
check_PROGRAMS += test_diagonalCovariance
check_PROGRAMS += test_fullCovariance
check_PROGRAMS += test_blockDiagonalCovariance
check_PROGRAMS += test_GslBlockMatrixInvertMultiply
check_PROGRAMS += test_unifiedPositionsOfMaximum
check_PROGRAMS += test_fullCovarianceRandomCoefficient
check_PROGRAMS += test_blockDiagonalCovarianceRandomCoefficients
check_PROGRAMS += test_fullCovarianceChain
check_PROGRAMS += test_diagonalCovarianceChain
check_PROGRAMS += test_scalarCovarianceChain
check_PROGRAMS += test_blockDiagonalCovarianceChain
check_PROGRAMS += test_1D_LinearLagrangeInterpolationSurrogate
check_PROGRAMS += test_2D_LinearLagrangeInterpolationSurrogate
check_PROGRAMS += test_3D_LinearLagrangeInterpolationSurrogate
check_PROGRAMS += test_4D_LinearLagrangeInterpolationSurrogate
check_PROGRAMS += test_InterpolationSurrogateHelper
check_PROGRAMS += test_build_InterpolationSurrogateBuilder
check_PROGRAMS += test_BoostInputOptionsParser
check_PROGRAMS += test_NoInputFile
check_PROGRAMS += test_optimizer_options
check_PROGRAMS += test_SharedPtr
check_PROGRAMS += test_serialEnv
check_PROGRAMS += test_jeffreys_pdf
check_PROGRAMS += test_gaussian_pdf_gradient
check_PROGRAMS += test_log_normal_pdf_gradient
check_PROGRAMS += test_beta_pdf_gradient
check_PROGRAMS += test_intercomm0_gravity
check_PROGRAMS += test_seq_of_vec_hdf5_write
check_PROGRAMS += test_optimizer_input_parameters
check_PROGRAMS += test_sip_gslopt_options

LDADD       = $(top_builddir)/src/libqueso.la

AM_CPPFLAGS = $(QUESO_CPPFLAGS)
AM_CPPFLAGS += -I.
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS +=  $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(ANN_CFLAGS)

if GRVY_ENABLED
  AM_CPPFLAGS += $(GRVY_CFLAGS)
endif

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE)
  LIBS += -lteuchoscore -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosremainder -lepetra
endif

if GLPK_ENABLED
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED
  AM_CPPFLAGS += $(HDF5_CFLAGS)
endif

if LIBMESH_ENABLED
  AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
endif

test_boxsubset_centroid_SOURCES = test_centroids/test_boxsubset_centroid.C
test_concatenation_centroid_SOURCES = test_centroids/test_concatenation_centroid.C
test_beta_mean_SOURCES = test_pdfmeans/test_beta_mean.C
test_concatenated_mean_SOURCES = test_pdfmeans/test_concatenated_mean.C
test_gamma_mean_SOURCES = test_pdfmeans/test_gamma_mean.C
test_uqEnvironmentNonFatal_SOURCES = test_Environment/test_uqEnvironmentNonFatal.C
test_uqGslVector_SOURCES = test_GslVector/test_uqGslVector.C
test_uqGaussianVectorRVClass_SOURCES = test_GaussianVectorRVClass/test_uqGaussianVectorRVClass.C
test_VectorPdf_gsl_SOURCES = test_GaussianVectorRVClass/test_VectorPdf_gsl.C
test_VectorRV_gsl_SOURCES = test_GaussianVectorRVClass/test_VectorRV_gsl.C
test_VectorRealizer_gsl_SOURCES = test_GaussianVectorRVClass/test_VectorRealizer_gsl.C
test_uqGslMatrix_SOURCES = test_GslMatrix/test_uqGslMatrix.C
test_uqTeuchosVector_SOURCES = test_TeuchosVector/test_uqTeuchosVector.C
test_uqexception_SOURCES = test_exception/test_exception.C
test_DistArrayCopy_SOURCES = test_DistArray/test_DistArrayCopy.C
test_DistArrayEquals_SOURCES = test_DistArray/test_DistArrayEquals.C
test_DistArrayMisc_SOURCES = test_DistArray/test_DistArrayMisc.C
test_GslVectorSpaceMisc_SOURCES = test_GslVectorSpace/test_GslVectorSpaceMisc.C
test_IntersectionSubsetContains_SOURCES = test_IntersectionSubset/test_IntersectionSubsetContains.C
test_operator_SOURCES = test_infinite/test_operator.C
test_inf_gaussian_SOURCES = test_infinite/test_inf_gaussian.C
test_inf_options_SOURCES = test_infinite/test_inf_options.C
test_SequenceOfVectorsErase_SOURCES = test_SequenceOfVectors/test_SequenceOfVectorsErase.C
test_GaussianMean1DRegression_SOURCES = test_Regression/test_GaussianMean1DRegression.C
test_gpmsa_cobra_SOURCES = test_Regression/test_gpmsa_cobra.C
test_gpmsa_vector_SOURCES = test_gpmsa/test_gpmsa_vector.C
test_LlhdTargetOutput_SOURCES = test_StatisticalInverseProblem/test_LlhdTargetOutput.C
test_jeffreys_SOURCES = test_Regression/test_jeffreys.C
test_gsloptimizer_SOURCES = test_optimizer/test_gsloptimizer.C
test_seedwithmap_SOURCES = test_optimizer/test_seedwithmap.C
test_seedwithmap_fd_SOURCES = test_optimizer/test_seedwithmap_fd.C
test_logitadaptedcov_SOURCES = test_Regression/test_logitadaptedcov.C
test_scalarCovariance_SOURCES = test_gaussian_likelihoods/test_scalarCovariance.C
test_diagonalCovariance_SOURCES = test_gaussian_likelihoods/test_diagonalCovariance.C
test_fullCovariance_SOURCES = test_gaussian_likelihoods/test_fullCovariance.C
test_blockDiagonalCovariance_SOURCES = test_gaussian_likelihoods/test_blockDiagonalCovariance.C
test_GslBlockMatrixInvertMultiply_SOURCES = test_GslBlockMatrix/test_GslBlockMatrixInvertMultiply.C
test_unifiedPositionsOfMaximum_SOURCES = test_SequenceOfVectors/test_unifiedPositionsOfMaximum.C
test_fullCovarianceRandomCoefficient_SOURCES = test_gaussian_likelihoods/test_fullCovarianceRandomCoefficient.C
test_blockDiagonalCovarianceRandomCoefficients_SOURCES = test_gaussian_likelihoods/test_blockDiagonalCovarianceRandomCoefficients.C
test_fullCovarianceChain_SOURCES = test_gaussian_likelihoods/test_fullCovarianceChain.C
test_diagonalCovarianceChain_SOURCES = test_gaussian_likelihoods/test_diagonalCovarianceChain.C
test_scalarCovarianceChain_SOURCES = test_gaussian_likelihoods/test_scalarCovarianceChain.C
test_blockDiagonalCovarianceChain_SOURCES = test_gaussian_likelihoods/test_blockDiagonalCovarianceChain.C
test_1D_LinearLagrangeInterpolationSurrogate_SOURCES = test_InterpolationSurrogate/test_1D_LinearLagrangeInterpolationSurrogate.C
test_2D_LinearLagrangeInterpolationSurrogate_SOURCES = test_InterpolationSurrogate/test_2D_LinearLagrangeInterpolationSurrogate.C
test_3D_LinearLagrangeInterpolationSurrogate_SOURCES = test_InterpolationSurrogate/test_3D_LinearLagrangeInterpolationSurrogate.C
test_4D_LinearLagrangeInterpolationSurrogate_SOURCES = test_InterpolationSurrogate/test_4D_LinearLagrangeInterpolationSurrogate.C
test_InterpolationSurrogateHelper_SOURCES = test_InterpolationSurrogate/test_InterpolationSurrogateHelper.C
test_build_InterpolationSurrogateBuilder_SOURCES = test_InterpolationSurrogate/test_build_InterpolationSurrogateBuilder.C
test_BoostInputOptionsParser_SOURCES = test_InputOptionsParser/test_BoostInputOptionsParser.C
test_NoInputFile_SOURCES = test_StatisticalInverseProblem/test_NoInputFile.C
test_optimizer_options_SOURCES = test_optimizer/test_optimizer_options.C
test_SharedPtr_SOURCES = pointers/test_SharedPtr.C
test_serialEnv_SOURCES = test_Environment/test_serialEnv.C
test_jeffreys_pdf_SOURCES = test_pdfs/test_jeffreys_pdf.C
test_gaussian_pdf_gradient_SOURCES = test_pdfs/test_gaussian_pdf_gradient.C
test_log_normal_pdf_gradient_SOURCES = test_pdfs/test_log_normal_pdf_gradient.C
test_beta_pdf_gradient_SOURCES = test_pdfs/test_beta_pdf_gradient.C

test_intercomm0_gravity_SOURCES =
test_intercomm0_gravity_SOURCES += test_intercomm0/gravity_likelihood.C
test_intercomm0_gravity_SOURCES += test_intercomm0/gravity_likelihood.h
test_intercomm0_gravity_SOURCES += test_intercomm0/gravity_main.C
test_intercomm0_gravity_SOURCES += test_intercomm0/gravity_qoi.C
test_intercomm0_gravity_SOURCES += test_intercomm0/gravity_qoi.h
test_intercomm0_gravity_CPPFLAGS = -I$(srcdir)/test_intercomm0 $(AM_CPPFLAGS)

test_seq_of_vec_hdf5_write_SOURCES = test_SequenceOfVectors/test_HDF5Write.C
test_optimizer_input_parameters_SOURCES = test_optimizer/test_optimizer_input_parameters.C
test_sip_gslopt_options_SOURCES = test_optimizer/test_sip_gslopt_options.C

# Files to freedom stamp
srcstamp =
srcstamp += $(test_boxsubset_centroid_SOURCES)
srcstamp += $(test_concatenation_centroid_SOURCES)
srcstamp += $(test_beta_mean_SOURCES)
srcstamp += $(test_concatenated_mean_SOURCES)
srcstamp += $(test_gamma_mean_SOURCES)
srcstamp += $(test_uqEnvironmentNonFatal_SOURCES)
srcstamp += $(test_uqEnvironmentOptionsPrint_SOURCES)
srcstamp += $(test_uqGslVectorConstructorFatal_SOURCES)
srcstamp += $(test_uqGslVector_SOURCES)
srcstamp += $(test_uqGaussianVectorRVClass_SOURCES)
srcstamp += $(test_VectorPdf_gsl_SOURCES)
srcstamp += $(test_VectorRV_gsl_SOURCES)
srcstamp += $(test_VectorRealizer_gsl_SOURCES)
srcstamp += $(test_uqGslMatrixConstructorFatal_SOURCES)
srcstamp += $(test_uqGslMatrix_SOURCES)
srcstamp += $(test_uqTeuchosVector_SOURCES)
srcstamp += $(test_uqexception_SOURCES)
srcstamp += $(test_DistArrayCtor_SOURCES)
srcstamp += $(test_DistArrayCopy_SOURCES)
srcstamp += $(test_DistArrayEquals_SOURCES)
srcstamp += $(test_DistArrayMisc_SOURCES)
srcstamp += $(test_GslVectorSpaceMisc_SOURCES)
srcstamp += $(test_IntersectionSubsetContains_SOURCES)
srcstamp += $(test_operator_SOURCES)
srcstamp += $(test_inf_gaussian_SOURCES)
srcstamp += $(test_inf_options_SOURCES)
srcstamp += $(test_GaussianMean1DRegression_SOURCES)
srcstamp += $(test_gpmsa_cobra_SOURCES)
srcstamp += $(test_gpmsa_vector_SOURCES)
srcstamp += $(test_LlhdTargetOutput_SOURCES)
srcstamp += $(test_jeffreys_SOURCES)
srcstamp += $(test_gsloptimizer_SOURCES)
srcstamp += $(test_seedwithmap_SOURCES)
srcstamp += $(test_seedwithmap_fd_SOURCES)
srcstamp += $(test_logitadaptedcov_SOURCES)
srcstamp += $(test_scalarCovariance_SOURCES)
srcstamp += $(test_diagonalCovariance_SOURCES)
srcstamp += $(test_fullCovariance_SOURCES)
srcstamp += $(test_blockDiagonalCovariance_SOURCES)
srcstamp += $(test_GslBlockMatrixInvertMultiply_SOURCES)
srcstamp += $(test_unifiedPositionsOfMaximum_SOURCES)
srcstamp += $(test_fullCovarianceRandomCoefficient_SOURCES)
srcstamp += $(test_blockDiagonalCovarianceRandomCoefficients_SOURCES)
srcstamp += $(test_1D_LinearLagrangeInterpolationSurrogate_SOURCES)
srcstamp += $(test_2D_LinearLagrangeInterpolationSurrogate_SOURCES)
srcstamp += $(test_3D_LinearLagrangeInterpolationSurrogate_SOURCES)
srcstamp += $(test_4D_LinearLagrangeInterpolationSurrogate_SOURCES)
srcstamp += $(test_InterpolationSurrogateHelper_SOURCES)
srcstamp += $(test_build_InterpolationSurrogateBuilder_SOURCES)
srcstamp += $(test_BoostInputOptionsParser_SOURCES)
srcstamp += $(test_NoInputFile_SOURCES)
srcstamp += $(test_optimizer_options_SOURCES)
srcstamp += $(test_SharedPtr_SOURCES)
srcstamp += $(test_serialEnv_SOURCES)
srcstamp += $(test_jeffreys_pdf_SOURCES)
srcstamp += $(test_gaussian_pdf_gradient_SOURCES)
srcstamp += $(test_log_normal_pdf_gradient_SOURCES)
srcstamp += $(test_beta_pdf_gradient_SOURCES)
srcstamp += $(test_intercomm0_gravity_SOURCES)
srcstamp += $(test_seq_of_vec_hdf5_write_SOURCES)
srcstamp += $(test_optimizer_input_parameters_SOURCES)
srcstamp += $(test_sip_gslopt_options_SOURCES)

TESTS =
TESTS += test_boxsubset_centroid
TESTS += test_concatenation_centroid
TESTS += test_beta_mean
TESTS += test_concatenated_mean
TESTS += test_gamma_mean
TESTS += test_uqEnvironmentNonFatal
TESTS += test_uqGslVector
TESTS += test_uqGaussianVectorRVClass
TESTS += test_VectorPdf_gsl
TESTS += test_VectorRV_gsl
TESTS += test_VectorRealizer_gsl
TESTS += test_uqGslMatrix
TESTS += test_uqTeuchosVector
TESTS += test_uqexception
TESTS += test_DistArrayCopy
TESTS += test_DistArrayEquals
TESTS += test_DistArrayMisc
TESTS += test_GslVectorSpaceMisc
TESTS += test_IntersectionSubsetContains
TESTS += test_operator
TESTS += test_inf_gaussian
TESTS += test_inf_options
TESTS += test_SequenceOfVectorsErase
TESTS += test_GaussianMean1DRegression
TESTS += test_Regression/test_cobra_samples_diff.sh
TESTS += test_gpmsa/test_gpmsa_samples_diff.sh
TESTS += test_StatisticalInverseProblem/test_LlhdTargetOutput.sh
TESTS += test_Regression/test_jeffreys_samples_diff.sh
TESTS += test_gsloptimizer
TESTS += test_seedwithmap
TESTS += test_seedwithmap_fd
TESTS += test_logitadaptedcov
TESTS += test_scalarCovariance
TESTS += test_diagonalCovariance
TESTS += test_fullCovariance
TESTS += test_blockDiagonalCovariance
TESTS += test_GslBlockMatrixInvertMultiply
TESTS += test_SequenceOfVectors/test_unifiedPositionsOfMaximum.sh
TESTS += test_fullCovarianceRandomCoefficient
TESTS += test_blockDiagonalCovarianceRandomCoefficients
TESTS += test_fullCovarianceChain
TESTS += test_diagonalCovarianceChain
TESTS += test_scalarCovarianceChain
TESTS += test_blockDiagonalCovarianceChain
TESTS += test_1D_LinearLagrangeInterpolationSurrogate
TESTS += test_2D_LinearLagrangeInterpolationSurrogate
TESTS += test_3D_LinearLagrangeInterpolationSurrogate
TESTS += test_4D_LinearLagrangeInterpolationSurrogate
TESTS += test_InterpolationSurrogateHelper
TESTS += test_build_InterpolationSurrogateBuilder
TESTS += test_BoostInputOptionsParser
TESTS += test_NoInputFile
TESTS += test_optimizer_options
TESTS += test_SharedPtr
TESTS += test_serialEnv
TESTS += test_jeffreys_pdf
TESTS += test_gaussian_pdf_gradient
TESTS += test_log_normal_pdf_gradient
TESTS += test_beta_pdf_gradient
TESTS += test_intercomm0/test_intercomm0_gravity_run.sh
TESTS += test_optimizer_input_parameters
TESTS += test_sip_gslopt_options
TESTS += test_SequenceOfVectors/test_seq_of_vec_hdf5_write_run.sh

XFAIL_TESTS = test_SequenceOfVectorsErase
if ! MPI_ENABLED
XFAIL_TESTS += test_SequenceOfVectors/test_unifiedPositionsOfMaximum.sh
endif


EXTRA_DIST =
EXTRA_DIST += common/compare.pl
EXTRA_DIST += common/verify.sh
EXTRA_DIST += test_infinite/inf_options
EXTRA_DIST += test_Regression/GaussianMean1DRegression_options
EXTRA_DIST += test_Regression/gpmsa_cobra_input.txt
EXTRA_DIST += test_Regression/test_cobra_samples_diff.sh
EXTRA_DIST += test_Regression/test_gpmsa_cobra_samples.m
EXTRA_DIST += test_Regression/ctf_dat.txt
EXTRA_DIST += test_Regression/dakota_pstudy.dat
EXTRA_DIST += test_gpmsa/gpmsa_vector_input.txt
EXTRA_DIST += test_gpmsa/test_gpmsa_samples_diff.sh
EXTRA_DIST += test_gpmsa/test_gpmsa_vector_samples.m
EXTRA_DIST += test_StatisticalInverseProblem/both_input.txt
EXTRA_DIST += test_StatisticalInverseProblem/llhdout_input.txt
EXTRA_DIST += test_StatisticalInverseProblem/targetout_input.txt
EXTRA_DIST += test_StatisticalInverseProblem/neither_input.txt
EXTRA_DIST += test_StatisticalInverseProblem/test_LlhdTargetOutput.sh
EXTRA_DIST += test_Regression/jeffreys_input.txt
EXTRA_DIST += test_Regression/test_jeffreys_samples_diff.sh
EXTRA_DIST += test_Regression/test_jeffreys_samples.m
EXTRA_DIST += test_Regression/adaptedcov_input.txt
EXTRA_DIST += test_gaussian_likelihoods/gaussian_consistency_input.txt
EXTRA_DIST += test_gaussian_likelihoods/queso_input.txt
EXTRA_DIST += test_InterpolationSurrogate/queso_input.txt
EXTRA_DIST += test_SequenceOfVectors/test_unifiedPositionsOfMaximum.sh
EXTRA_DIST += test_InputOptionsParser/test_options_good.txt
EXTRA_DIST += test_InputOptionsParser/test_options_bad.txt
EXTRA_DIST += test_InputOptionsParser/test_options_default.txt
EXTRA_DIST += test_optimizer/input_test_optimizer_options
EXTRA_DIST += test_Environment/input_test_serialEnv
EXTRA_DIST += test_intercomm0/gravity_1proc.txt
EXTRA_DIST += test_intercomm0/gravity_2proc.txt
EXTRA_DIST += test_intercomm0/test_intercomm0_gravity_run.sh
EXTRA_DIST += test_optimizer/input_test_optimizer_input_parameters
EXTRA_DIST += test_SequenceOfVectors/test_seq_of_vec_hdf5_write_run.sh

CLEANFILES =
CLEANFILES += test_Environment/debug_output_sub0.txt
CLEANFILES += gslvector_out_sub0.m
CLEANFILES += test_write_InterpolationSurrogateBuilder_1.dat
CLEANFILES += test_write_InterpolationSurrogateBuilder_2.dat

clean-local:
	rm -rf $(top_builddir)/test/chain0
	rm -rf $(top_builddir)/test/outputData
	rm -rf $(top_builddir)/test/test_output_gaussian_likelihoods
	rm -rf $(top_builddir)/test/test_gpmsa_cobra_output
	rm -rf $(top_builddir)/test/test_gpmsa_vector_output
	rm -rf $(top_builddir)/test/test_adaptedcov_output
	rm -rf $(top_builddir)/test/test_outputNoInputFile
	rm -rf $(top_builddir)/test/test_output_interp_surrogates
	rm -rf $(top_builddir)/test/output_test_optimizer_options
	rm -rf $(top_builddir)/test/output_test_serialEnv
	rm -rf $(top_builddir)/test/output_test_intercomm0_gravity_1
	rm -rf $(top_builddir)/test/output_test_intercomm0_gravity_2

if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

if GIT_CLONE
.license.stamp: $(top_srcdir)/LICENSE
	@$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(QUESO_STAMPED_FILES)
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
