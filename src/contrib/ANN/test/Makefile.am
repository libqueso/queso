## Process this file with automake to produce Makefile.in
if HAVE_ANN

AM_CPPFLAGS = $(ANN_CFLAGS)

noinst_PROGRAMS = ann_test

ann_test_SOURCES =
ann_test_SOURCES += rand.cpp
ann_test_SOURCES += ann_test.cpp
ann_test_LDADD = $(top_builddir)/src/contrib/ANN/libANN.la

# Required for AX_AM_MACROS
@INC_AMINCLUDE@

endif
