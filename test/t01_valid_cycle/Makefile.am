if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

AM_CPPFLAGS  = -I. \
	-I$(top_srcdir)/src/core/inc \
	-I$(top_srcdir)/src/misc/inc \
	-I$(top_srcdir)/src/basic/inc \
	-I$(top_srcdir)/src/stats/inc \
	-I$(TRILINOS_INCLUDE) \
          $(EPETRA_INCLUDES) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(GRVY_CFLAGS) $(GLPK_CFLAGS) $(HDF5_CFLAGS)

LIBS   = -L$(top_builddir)/src/ -lqueso

EXTRA_DIST     = inputData regression rtest01.sh
TESTS          = rtest01.sh 
check_PROGRAMS = TgaValidationCycle_gsl

TgaValidationCycle_gsl_SOURCES = TgaValidationCycle_gsl.C TgaValidationCycle_appl.h \
				TgaValidationCycle_likelihood.h TgaValidationCycle_qoi.h \
				tgaCycle.inp

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.C *.h
	echo 'updated source license headers' >$@
endif

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

# Required for AX_AM_MACROS
@INC_AMINCLUDE@



