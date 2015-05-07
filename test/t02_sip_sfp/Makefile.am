if GIT_CLONE
BUILT_SOURCES = .license.stamp
endif

AM_CPPFLAGS =
AM_CPPFLAGS += -I.
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS +=  $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(ANN_CFLAGS)

LIBS   = -L$(top_builddir)/src/ -lqueso

if GRVY_ENABLED
  AM_CPPFLAGS += $(GRVY_CFLAGS)
endif

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE)
  LIBS += -lteuchoscore -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosremainder -lepetra
endif

if GLPK_ENABLED
  AM_CPPFLAGS += $(GLPK_CFLAGS)
  LIBS += $(GLPK_LIBS)
endif

if HDF5_ENABLED
  AM_CPPFLAGS += $(HDF5_CFLAGS)
  LIBS += $(HDF5_LIBS)
endif

EXTRA_DIST     = regression rtest02.sh
check_PROGRAMS = SipSfpExample_gsl
TESTS          = rtest02.sh

SipSfpExample_gsl_SOURCES =
SipSfpExample_gsl_SOURCES += example_main.C example_compute.C example_likelihood.C example_qoi.C
SipSfpExample_gsl_SOURCES += example_compute.h example_likelihood.h example_qoi.h
SipSfpExample_gsl_SOURCES += example.inp

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

if GIT_CLONE
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(QUESO_STAMPED_FILES)
	echo 'updated source license headers' >$@
endif

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
