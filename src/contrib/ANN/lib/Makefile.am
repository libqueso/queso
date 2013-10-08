## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES		= libANN.la
libANN_includedir 	= $(prefix)/include/ANN

ANN_DIR = $(top_srcdir)/src/contrib/ANN

AM_CPPFLAGS = -I$(ANN_DIR)/src -I$(ANN_DIR)/include

libANN_la_SOURCES =
libANN_la_SOURCES += $(ANN_DIR)/src/ANN.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/brute.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_tree.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_util.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_split.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_dump.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_search.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_pr_search.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/kd_fix_rad_search.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/bd_tree.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/bd_search.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/bd_pr_search.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/bd_fix_rad_search.cpp
libANN_la_SOURCES += $(ANN_DIR)/src/perf.cpp

libANN_include_HEADERS =
libANN_include_HEADERS += $(ANN_DIR)/include/ANN/ANN.h
libANN_include_HEADERS += $(ANN_DIR)/include/ANN/ANNx.h
libANN_include_HEADERS += $(ANN_DIR)/include/ANN/ANNperf.h
