SUBDIRS    = t01_valid_cycle t02_sip_sfp t03_sequence t04_bimodal

if GIT_CLONE
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += test_uqEnvironmentCopy
check_PROGRAMS += test_uqEnvironmentEquals
check_PROGRAMS += test_uqEnvironmentNonFatal
check_PROGRAMS += test_uqEnvironmentOptionsPrint
check_PROGRAMS += test_uqGslVectorConstructorFatal
check_PROGRAMS += test_uqGslVector
check_PROGRAMS += test_uqGaussianVectorRVClass
check_PROGRAMS += test_VectorPdf_gsl
check_PROGRAMS += test_VectorRV_gsl
check_PROGRAMS += test_VectorRealizer_gsl
check_PROGRAMS += test_uqGslMatrixConstructorFatal
check_PROGRAMS += test_uqGslMatrix
check_PROGRAMS += test_uqTeuchosVector
check_PROGRAMS += test_uqexception
check_PROGRAMS += test_DistArrayCtor
check_PROGRAMS += test_DistArrayCopy
check_PROGRAMS += test_DistArrayEquals
check_PROGRAMS += test_DistArrayMisc
check_PROGRAMS += test_GslVectorSpaceCtor
check_PROGRAMS += test_GslVectorSpaceMisc
check_PROGRAMS += test_IntersectionSubsetContains
check_PROGRAMS += test_operator
check_PROGRAMS += test_inf_gaussian
check_PROGRAMS += test_inf_options
check_PROGRAMS += test_SequenceOfVectorsErase
check_PROGRAMS += test_GaussianMean1DRegression
check_PROGRAMS += test_gpmsa_cobra
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

LDADD       = $(top_builddir)/src/libqueso.la

AM_CPPFLAGS =
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

test_uqEnvironmentCopy_SOURCES = test_Environment/test_uqEnvironmentCopy.C
test_uqEnvironmentEquals_SOURCES = test_Environment/test_uqEnvironmentEquals.C
test_uqEnvironmentNonFatal_SOURCES = test_Environment/test_uqEnvironmentNonFatal.C
test_uqEnvironmentOptionsPrint_SOURCES = test_uqEnvironmentOptions/test_uqEnvironmentOptionsPrint.C
test_uqGslVectorConstructorFatal_SOURCES = test_GslVector/test_uqGslVectorConstructorFatal.C
test_uqGslVector_SOURCES = test_GslVector/test_uqGslVector.C
test_uqGaussianVectorRVClass_SOURCES = test_GaussianVectorRVClass/test_uqGaussianVectorRVClass.C
test_VectorPdf_gsl_SOURCES = test_GaussianVectorRVClass/test_VectorPdf_gsl.C
test_VectorRV_gsl_SOURCES = test_GaussianVectorRVClass/test_VectorRV_gsl.C
test_VectorRealizer_gsl_SOURCES = test_GaussianVectorRVClass/test_VectorRealizer_gsl.C
test_uqGslMatrixConstructorFatal_SOURCES = test_GslMatrix/test_uqGslMatrixConstructorFatal.C
test_uqGslMatrix_SOURCES = test_GslMatrix/test_uqGslMatrix.C
test_uqTeuchosVector_SOURCES = test_TeuchosVector/test_uqTeuchosVector.C
test_uqexception_SOURCES = test_exception/test_exception.C
test_DistArrayCtor_SOURCES = test_DistArray/test_DistArrayCtor.C
test_DistArrayCopy_SOURCES = test_DistArray/test_DistArrayCopy.C
test_DistArrayEquals_SOURCES = test_DistArray/test_DistArrayEquals.C
test_DistArrayMisc_SOURCES = test_DistArray/test_DistArrayMisc.C
test_GslVectorSpaceCtor_SOURCES = test_GslVectorSpace/test_GslVectorSpaceCtor.C
test_GslVectorSpaceMisc_SOURCES = test_GslVectorSpace/test_GslVectorSpaceMisc.C
test_IntersectionSubsetContains_SOURCES = test_IntersectionSubset/test_IntersectionSubsetContains.C
test_operator_SOURCES = test_infinite/test_operator.C
test_inf_gaussian_SOURCES = test_infinite/test_inf_gaussian.C
test_inf_options_SOURCES = test_infinite/test_inf_options.C
test_SequenceOfVectorsErase_SOURCES = test_SequenceOfVectors/test_SequenceOfVectorsErase.C
test_GaussianMean1DRegression_SOURCES = test_Regression/test_GaussianMean1DRegression.C
test_GaussianMean1DRegression_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
test_gpmsa_cobra_SOURCES = test_Regression/test_gpmsa_cobra.C
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

# Files to freedom stamp
srcstamp =
srcstamp += $(test_uqEnvironmentCopy_SOURCES)
srcstamp += $(test_uqEnvironmentEquals_SOURCES)
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
srcstamp += $(test_GslVectorSpaceCtor_SOURCES)
srcstamp += $(test_GslVectorSpaceMisc_SOURCES)
srcstamp += $(test_IntersectionSubsetContains_SOURCES)
srcstamp += $(test_operator_SOURCES)
srcstamp += $(test_inf_gaussian_SOURCES)
srcstamp += $(test_inf_options_SOURCES)
srcstamp += $(test_GaussianMean1DRegression_SOURCES)
srcstamp += $(test_gpmsa_cobra_SOURCES)
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

TESTS =
TESTS += $(top_builddir)/test/test_uqEnvironmentCopy
TESTS += $(top_builddir)/test/test_uqEnvironmentEquals
TESTS += $(top_builddir)/test/test_uqEnvironmentNonFatal
TESTS += $(top_builddir)/test/test_uqEnvironmentOptionsPrint
TESTS += $(top_builddir)/test/test_uqGslVectorConstructorFatal
TESTS += $(top_builddir)/test/test_uqGslVector
TESTS += $(top_builddir)/test/test_uqGaussianVectorRVClass
TESTS += $(top_builddir)/test/test_VectorPdf_gsl
TESTS += $(top_builddir)/test/test_VectorRV_gsl
TESTS += $(top_builddir)/test/test_VectorRealizer_gsl
TESTS += $(top_builddir)/test/test_uqGslMatrixConstructorFatal
TESTS += $(top_builddir)/test/test_uqGslMatrix
TESTS += $(top_builddir)/test/test_uqTeuchosVector
TESTS += $(top_builddir)/test/test_uqexception
TESTS += $(top_builddir)/test/test_DistArrayCtor
TESTS += $(top_builddir)/test/test_DistArrayCopy
TESTS += $(top_builddir)/test/test_DistArrayEquals
TESTS += $(top_builddir)/test/test_DistArrayMisc
TESTS += $(top_builddir)/test/test_GslVectorSpaceCtor
TESTS += $(top_builddir)/test/test_GslVectorSpaceMisc
TESTS += $(top_builddir)/test/test_IntersectionSubsetContains
TESTS += $(top_builddir)/test/test_operator
TESTS += $(top_builddir)/test/test_inf_gaussian
TESTS += $(top_builddir)/test/test_inf_options
TESTS += $(top_builddir)/test/test_SequenceOfVectorsErase
TESTS += $(top_builddir)/test/test_GaussianMean1DRegression
TESTS += $(top_builddir)/test/test_Regression/test_cobra_samples_diff.sh
TESTS += $(top_builddir)/test/test_StatisticalInverseProblem/test_LlhdTargetOutput.sh
TESTS += $(top_builddir)/test/test_Regression/test_jeffreys_samples_diff.sh
TESTS += $(top_builddir)/test/test_gsloptimizer
TESTS += $(top_builddir)/test/test_seedwithmap
TESTS += $(top_builddir)/test/test_seedwithmap_fd
TESTS += $(top_builddir)/test/test_logitadaptedcov
TESTS += $(top_builddir)/test/test_scalarCovariance
TESTS += $(top_builddir)/test/test_diagonalCovariance
TESTS += $(top_builddir)/test/test_fullCovariance
TESTS += $(top_builddir)/test/test_blockDiagonalCovariance
TESTS += $(top_builddir)/test/test_GslBlockMatrixInvertMultiply
TESTS += $(top_builddir)/test/test_SequenceOfVectors/test_unifiedPositionsOfMaximum.sh

XFAIL_TESTS = $(top_builddir)/test/test_SequenceOfVectorsErase

EXTRA_DIST =
EXTRA_DIST += common/compare.pl
EXTRA_DIST += common/verify.sh
EXTRA_DIST += test_uqEnvironmentOptions/test.inp
EXTRA_DIST += test_Environment/copy_env
EXTRA_DIST += test_infinite/inf_options
EXTRA_DIST += test_Regression/GaussianMean1DRegression_options.in
EXTRA_DIST += test_Regression/gpmsa_cobra_input.txt.in
EXTRA_DIST += test_Regression/test_cobra_samples_diff.sh.in
EXTRA_DIST += test_Regression/test_gpmsa_cobra_samples.m.in
EXTRA_DIST += test_Regression/ctf_dat.txt.in
EXTRA_DIST += test_Regression/dakota_pstudy.dat.in
EXTRA_DIST += test_StatisticalInverseProblem/both_input.txt.in
EXTRA_DIST += test_StatisticalInverseProblem/llhdout_input.txt.in
EXTRA_DIST += test_StatisticalInverseProblem/targetout_input.txt.in
EXTRA_DIST += test_StatisticalInverseProblem/neither_input.txt.in
EXTRA_DIST += test_StatisticalInverseProblem/test_LlhdTargetOutput.sh.in
EXTRA_DIST += test_Regression/jeffreys_input.txt.in
EXTRA_DIST += test_Regression/test_jeffreys_samples_diff.sh.in
EXTRA_DIST += test_Regression/test_jeffreys_samples.m.in
EXTRA_DIST += test_Regression/adaptedcov_input.txt.in
EXTRA_DIST += test_gaussian_likelihoods/queso_input.txt.in
EXTRA_DIST += test_SequenceOfVectors/test_unifiedPositionsOfMaximum.sh.in

CLEANFILES =
CLEANFILES += test_Environment/debug_output_sub0.txt
CLEANFILES += gslvector_out_sub0.m

clean-local:
	rm -rf $(top_builddir)/test/chain0
	rm -rf $(top_builddir)/test/outputData
	rm -rf $(top_builddir)/test/test_output_gaussian_likelihoods
	rm -rf $(top_builddir)/test/test_gpmsa_cobra_output

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
