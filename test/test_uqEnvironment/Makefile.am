if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += test_uqEnvironment
check_PROGRAMS += test_uqEnvironmentCopy
check_PROGRAMS += test_uqEnvironmentEquals
check_PROGRAMS += test_uqEnvironmentNonFatal

LIBS         = -L$(top_builddir)/src/ -lqueso 

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src \
              -I. \
              -I$(top_srcdir)/src/core/inc \
              -I$(top_srcdir)/src/misc/inc \
              -I$(top_srcdir)/src/basic/inc \
              -I$(top_srcdir)/src/stats/inc \
              $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(GRVY_CFLAGS) $(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
endif

test_uqEnvironment_SOURCES = $(top_srcdir)/test/test_uqEnvironment/test_uqEnvironment.C
test_uqEnvironmentCopy_SOURCES = $(top_srcdir)/test/test_uqEnvironment/test_uqEnvironmentCopy.C
test_uqEnvironmentEquals_SOURCES = $(top_srcdir)/test/test_uqEnvironment/test_uqEnvironmentEquals.C
test_uqEnvironmentNonFatal_SOURCES = $(top_srcdir)/test/test_uqEnvironment/test_uqEnvironmentNonFatal.C

TESTS = $(top_builddir)/test/test_uqEnvironment/test_uqEnvironment.sh \
				$(top_builddir)/test/test_uqEnvironment/test_uqEnvironmentCopy.sh \
				$(top_builddir)/test/test_uqEnvironment/test_uqEnvironmentEquals.sh \
				$(top_builddir)/test/test_uqEnvironment/test_uqEnvironmentNonFatal

EXTRA_DIST = test_uqEnvironment.sh \
             test_uqEnvironmentCopy.sh \
             test_uqEnvironmentEquals.sh \
						 test_uqEnvironmentNonFatal

CLEANFILES = $(top_srcdir)/test/test_uqEnvironment/debug_output_sub0.txt
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.C
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
