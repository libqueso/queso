if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

check_PROGRAMS = # Append below
check_PROGRAMS += test_uqEnvironment

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

test_uqEnvironment_SOURCES = test_uqEnvironment.C

TESTS = test_uqEnvironment.sh
EXTRA_DIST = test_uqEnvironment.sh

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.C
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
