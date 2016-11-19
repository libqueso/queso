BUILT_SOURCES =

AM_CPPFLAGS = $(QUESO_CPPFLAGS)
AM_CPPFLAGS += -I.
AM_CPPFLAGS += -I$(top_builddir)/inc
AM_CPPFLAGS += -I$(TRILINOS_INCLUDE)
AM_CPPFLAGS += $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS) $(GLPK_CFLAGS) $(HDF5_CFLAGS)

LIBS  = -L$(top_builddir)/src/ -lqueso

EXTRA_DIST                  = example_compute.h
noinst_PROGRAMS             = SequenceExample_gsl
SequenceExample_gsl_SOURCES = example_main.C example_compute.C

if QUESO_LICENSESTAMPEXISTS
BUILT_SOURCES += .license.stamp
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(QUESO_STAMPED_FILES)
	echo 'updated source license headers' >$@
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
