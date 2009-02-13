## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/libs/basic/inc $(MPI_CXXFLAGS) $(BOOST_CPPFLAGS)

lib_LIBRARIES = 

lib_LIBRARIES += libqueso_stats_any.a
libqueso_stats_any_a_SOURCES = 

if UQBT_GSL
lib_LIBRARIES += libqueso_stats_gsl.a
libqueso_stats_gsl_a_SOURCES =
endif