if SVN_CHECKOUT
BUILT_SOURCES = .license.stamp
endif

AM_CPPFLAGS  = -I. \
	-I$(top_srcdir)/src/core/inc \
	-I$(top_srcdir)/src/misc/inc \
	-I$(top_srcdir)/src/basic/inc \
	-I$(top_srcdir)/src/stats/inc \
         $(BOOST_CPPFLAGS) $(GSL_CFLAGS) $(GRVY_CFLAGS) $(ANN_CFLAGS)

if TRILINOS_ENABLED
  AM_CPPFLAGS += -I$(TRILINOS_INCLUDE) $(EPETRA_INCLUDES)
  LIBS += -lteuchos -lepetra
endif

if GLPK_ENABLED	
  AM_CPPFLAGS += $(GLPK_CFLAGS)
  LIBS += $(GLPK_LIBS)
endif

if HDF5_ENABLED	
  AM_CPPFLAGS += $(HDF5_CFLAGS) 
  LIBS += $(HDF5_LIBS)
endif

LIBS  = -L$(top_builddir)/src/ -lqueso

EXTRA_DIST     = regression rtest04.sh
check_PROGRAMS = BimodalExample_gsl
TESTS          = rtest04.sh 

BimodalExample_gsl_SOURCES = example_main.C example_compute.C example_likelihood.C \
			example_compute.h example_likelihood.h \
			example_1chain.inp

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

if SVN_CHECKOUT
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(srcdir)/*.C $(srcdir)/*.h
	echo 'updated source license headers' >$@
endif

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

