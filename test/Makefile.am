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

LIBS         = -L$(top_builddir)/src/ -lqueso

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

test_uqEnvironmentCopy_SOURCES = $(top_srcdir)/test/test_Environment/test_uqEnvironmentCopy.C
test_uqEnvironmentEquals_SOURCES = $(top_srcdir)/test/test_Environment/test_uqEnvironmentEquals.C
test_uqEnvironmentNonFatal_SOURCES = $(top_srcdir)/test/test_Environment/test_uqEnvironmentNonFatal.C
test_uqEnvironmentOptionsPrint_SOURCES = $(top_srcdir)/test/test_uqEnvironmentOptions/test_uqEnvironmentOptionsPrint.C
test_uqGslVectorConstructorFatal_SOURCES = $(top_srcdir)/test/test_GslVector/test_uqGslVectorConstructorFatal.C
test_uqGslVector_SOURCES = $(top_srcdir)/test/test_GslVector/test_uqGslVector.C
test_uqGaussianVectorRVClass_SOURCES = $(top_srcdir)/test/test_GaussianVectorRVClass/test_uqGaussianVectorRVClass.C
test_uqGslMatrixConstructorFatal_SOURCES = $(top_srcdir)/test/test_GslMatrix/test_uqGslMatrixConstructorFatal.C
test_uqGslMatrix_SOURCES = $(top_srcdir)/test/test_GslMatrix/test_uqGslMatrix.C
test_uqTeuchosVector_SOURCES = $(top_srcdir)/test/test_TeuchosVector/test_uqTeuchosVector.C
test_uqexception_SOURCES = $(top_srcdir)/test/test_exception/test_exception.C
test_DistArrayCtor_SOURCES = $(top_srcdir)/test/test_DistArray/test_DistArrayCtor.C
test_DistArrayCopy_SOURCES = $(top_srcdir)/test/test_DistArray/test_DistArrayCopy.C
test_DistArrayEquals_SOURCES = $(top_srcdir)/test/test_DistArray/test_DistArrayEquals.C
test_DistArrayMisc_SOURCES = $(top_srcdir)/test/test_DistArray/test_DistArrayMisc.C
test_GslVectorSpaceCtor_SOURCES = $(top_srcdir)/test/test_GslVectorSpace/test_GslVectorSpaceCtor.C
test_GslVectorSpaceMisc_SOURCES = $(top_srcdir)/test/test_GslVectorSpace/test_GslVectorSpaceMisc.C
test_IntersectionSubsetContains_SOURCES = $(top_srcdir)/test/test_IntersectionSubset/test_IntersectionSubsetContains.C
test_operator_SOURCES = $(top_srcdir)/test/test_infinite/test_operator.C
test_inf_gaussian_SOURCES = $(top_srcdir)/test/test_infinite/test_inf_gaussian.C
test_inf_options_SOURCES = $(top_srcdir)/test/test_infinite/test_inf_options.C
test_SequenceOfVectorsErase_SOURCES = $(top_srcdir)/test/test_SequenceOfVectors/test_SequenceOfVectorsErase.C
test_GaussianMean1DRegression_SOURCES = $(top_srcdir)/test/test_Regression/test_GaussianMean1DRegression.C

# Files to freedom stamp
srcstamp =
srcstamp += $(test_uqEnvironmentCopy_SOURCES)
srcstamp += $(test_uqEnvironmentEquals_SOURCES)
srcstamp += $(test_uqEnvironmentNonFatal_SOURCES)
srcstamp += $(test_uqEnvironmentOptionsPrint_SOURCES)
srcstamp += $(test_uqGslVectorConstructorFatal_SOURCES)
srcstamp += $(test_uqGslVector_SOURCES)
srcstamp += $(test_uqGaussianVectorRVClass_SOURCES)
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

TESTS =
TESTS += $(top_builddir)/test/test_uqEnvironmentCopy
TESTS += $(top_builddir)/test/test_uqEnvironmentEquals
TESTS += $(top_builddir)/test/test_uqEnvironmentNonFatal
TESTS += $(top_builddir)/test/test_uqEnvironmentOptionsPrint
TESTS += $(top_builddir)/test/test_uqGslVectorConstructorFatal
TESTS += $(top_builddir)/test/test_uqGslVector
TESTS += $(top_builddir)/test/test_uqGaussianVectorRVClass
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

XFAIL_TESTS = $(top_builddir)/test/test_SequenceOfVectorsErase

EXTRA_DIST =
EXTRA_DIST += common/compare.pl
EXTRA_DIST += common/verify.sh
EXTRA_DIST += test_uqEnvironmentOptions/test.inp
EXTRA_DIST += test_Environment/copy_env
EXTRA_DIST += test_infinite/inf_options
EXTRA_DIST += test_Regression/test_GaussianMean1DRegression.inp

CLEANFILES =
CLEANFILES += $(top_srcdir)/test/test_Environment/debug_output_sub0.txt
CLEANFILES += $(top_srcdir)/test/gslvector_out_sub0.m

clean-local:
	rm -rf $(top_builddir)/test/chain0
	rm -rf $(top_builddir)/test/outputData

if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

if GIT_CLONE
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(srcstamp)
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
