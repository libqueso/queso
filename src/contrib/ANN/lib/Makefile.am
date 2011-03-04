## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES		= libANN.la
libANN_includedir 	= $(prefix)/include/ANN

ANN_DIR = $(top_srcdir)/src/contrib/ANN

AM_CPPFLAGS = \
	-I$(ANN_DIR)/src       \
	-I$(ANN_DIR)/include

libANN_la_SOURCES = \
	$(ANN_DIR)/src/ANN.cpp \
	$(ANN_DIR)/src/brute.cpp \
	$(ANN_DIR)/src/kd_tree.cpp \
	$(ANN_DIR)/src/kd_util.cpp \
	$(ANN_DIR)/src/kd_split.cpp \
	$(ANN_DIR)/src/kd_dump.cpp \
	$(ANN_DIR)/src/kd_search.cpp \
	$(ANN_DIR)/src/kd_pr_search.cpp \
	$(ANN_DIR)/src/kd_fix_rad_search.cpp \
	$(ANN_DIR)/src/bd_tree.cpp \
	$(ANN_DIR)/src/bd_search.cpp \
	$(ANN_DIR)/src/bd_pr_search.cpp \
	$(ANN_DIR)/src/bd_fix_rad_search.cpp \
	$(ANN_DIR)/src/perf.cpp

libANN_include_HEADERS = \
	$(ANN_DIR)/include/ANN/ANN.h \
	$(ANN_DIR)/include/ANN/ANNx.h \
	$(ANN_DIR)/include/ANN/ANNperf.h
