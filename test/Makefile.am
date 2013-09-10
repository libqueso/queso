SUBDIRS    = t01_valid_cycle t02_sip_sfp t03_sequence t04_bimodal

if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += test_uqEnvironment
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

LIBS         = -L$(top_builddir)/src/ -lqueso

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src \
              -I. \
              -I$(top_srcdir)/src/core/inc \
              -I$(top_srcdir)/src/misc/inc \
              -I$(top_srcdir)/src/basic/inc \
              -I$(top_srcdir)/src/stats/inc \
              -I$(top_srcdir)/src/gp/inc \
              $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(ANN_CFLAGS)

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

test_uqEnvironment_SOURCES = $(top_srcdir)/test/test_Environment/test_uqEnvironment.C
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

# Files to freedom stamp
srcstamp = $(test_uqEnvironment_SOURCES) \
					 $(test_uqEnvironmentCopy_SOURCES) \
					 $(test_uqEnvironmentEquals_SOURCES) \
					 $(test_uqEnvironmentNonFatal_SOURCES) \
					 $(test_uqEnvironmentOptionsPrint_SOURCES) \
					 $(test_uqGslVectorConstructorFatal_SOURCES) \
					 $(test_uqGslVector_SOURCES) \
					 $(test_uqGaussianVectorRVClass_SOURCES) \
           $(test_uqGslMatrixConstructorFatal_SOURCES) \
					 $(test_uqGslMatrix_SOURCES)


TESTS = $(top_builddir)/test/test_Environment/test_uqEnvironment.sh \
				$(top_builddir)/test/test_Environment/test_uqEnvironmentCopy.sh \
				$(top_builddir)/test/test_Environment/test_uqEnvironmentEquals.sh \
				$(top_builddir)/test/test_uqEnvironmentNonFatal \
        $(top_builddir)/test/test_uqEnvironmentOptions/test_uqEnvironmentOptionsPrint.sh \
				$(top_builddir)/test/test_GslVector/test_uqGslVectorConstructorFatal.sh \
				$(top_builddir)/test/test_uqGslVector \
        $(top_builddir)/test/test_uqGaussianVectorRVClass \
				$(top_builddir)/test/test_GslMatrix/test_uqGslMatrixConstructorFatal.sh \
				$(top_builddir)/test/test_uqGslMatrix \
				$(top_builddir)/test/test_uqTeuchosVector
				$(top_builddir)/test/test_exception

EXTRA_DIST = common/compare.pl \
						 common/verify.sh \
             test_Environment/test_uqEnvironment.sh \
             test_Environment/test_uqEnvironmentCopy.sh \
             test_Environment/test_uqEnvironmentEquals.sh \
						 test_GslVector/test_uqGslVectorConstructorFatal.sh \
             test_GslMatrix/test_uqGslMatrixConstructorFatal.sh \
						 test_uqEnvironmentOptions/test_uqEnvironmentOptionsPrint.sh \
						 test_uqEnvironmentOptions/test.inp

CLEANFILES = $(top_srcdir)/test/test_Environment/debug_output_sub0.txt \
						 $(top_srcdir)/test/gslvector_out_sub0.m

if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(srcstamp)
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
