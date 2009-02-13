## Process this file with automake to produce Makefile.in

# Build in these directories:
DIST_SUBDIRS = basic stats

AM_CPPFLAGS = -I$(top_srcdir)/libs/basic/inc -I$(top_srcdir)/libs/stats/inc $(MPI_CXXFLAGS) $(BOOST_CPPFLAGS)

lib_LIBRARIES = libqueso.a

# Sources from basic/src
libqueso_a_SOURCES = $(top_srcdir)/libs/basic/src/uqEnvironment.C \
	$(top_srcdir)/libs/basic/src/uqDefines.C \
	$(top_srcdir)/libs/basic/src/uqVector.C \
	$(top_srcdir)/libs/basic/src/uqMatrix.C \
	$(top_srcdir)/libs/basic/src/uqMiscellaneous.C \
	$(top_srcdir)/libs/basic/src/uqRealFft.C \
	$(top_srcdir)/libs/basic/src/uqComplexFft.C \
	$(top_srcdir)/libs/basic/src/uqChainStatisticalOptions.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVector.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosMatrix.C \
	$(top_srcdir)/libs/basic/src/uqTrilinosVectorSpace.C

# Sources from basic/src added conditional on gsl
if UQBT_GSL
libqueso_a_SOURCES += $(top_srcdir)/libs/basic/src/uqGslVector.C \
	$(top_srcdir)/libs/basic/src/uqGslMatrix.C \
	$(top_srcdir)/libs/basic/src/uqGslOdeSolver.C \
	$(top_srcdir)/libs/basic/src/uqGslVectorSpace.C
endif

# Sources from stats/src
# Currently none... add when necessary