## Process this file with automake to produce Makefile.in
if HAVE_ANN

noinst_LTLIBRARIES = libANN.la

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include -I$(top_srcdir)/src/contrib/ANN/include

libANN_la_SOURCES =
libANN_la_SOURCES += src/ANN.cpp
libANN_la_SOURCES += src/brute.cpp
libANN_la_SOURCES += src/kd_tree.cpp
libANN_la_SOURCES += src/kd_util.cpp
libANN_la_SOURCES += src/kd_split.cpp
libANN_la_SOURCES += src/kd_dump.cpp
libANN_la_SOURCES += src/kd_search.cpp
libANN_la_SOURCES += src/kd_pr_search.cpp
libANN_la_SOURCES += src/kd_fix_rad_search.cpp
libANN_la_SOURCES += src/bd_tree.cpp
libANN_la_SOURCES += src/bd_search.cpp
libANN_la_SOURCES += src/bd_pr_search.cpp
libANN_la_SOURCES += src/bd_fix_rad_search.cpp
libANN_la_SOURCES += src/perf.cpp
libANN_la_SOURCES += include/ANN/ANN.h
libANN_la_SOURCES += include/ANN/ANNx.h
libANN_la_SOURCES += include/ANN/ANNperf.h

endif
