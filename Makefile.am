## Process this file with automake to produce Makefile.in

include $(top_srcdir)/doxygen/aminclude.am

# Don't require GNU-style project
AUTOMAKE_OPTIONS = foreign

# Build in these directories:
SUBDIRS = libs \
	examples \
	test/t01_valid_cycle/valid_cycle \
	test/t02_sip_sfp/sip_sfp \
	test/t03_sequence/sequence \
        test/gsl_tests \
        test

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/missing

EXTRA_DIST = $(DX_CONFIG) \
	queso.page

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES =     \
	aclocal.m4             \
	autom4te.cache/*       \
	$(AUX_DIST)            \
	config.log             \
	config.status          \
	config.sub             \
	configure              \
	config.h	       \
	config.h.in	       \
	Makefile.in            \
	libs/Makefile.in       \
	libs/core/Makefile.in \
	libs/core/src/Makefile.in \
	libs/core/inc/Makefile.in \
	libs/misc/Makefile.in \
	libs/misc/src/Makefile.in \
	libs/misc/inc/Makefile.in \
	libs/basic/Makefile.in \
	libs/basic/src/Makefile.in \
	libs/basic/inc/Makefile.in \
	libs/stats/Makefile.in \
	libs/stats/src/Makefile.in \
	libs/stats/inc/Makefile.in \
	examples/Makefile.in \
	examples/markovChain1/Makefile.in \
	examples/markovChain1/src/Makefile.in \
	examples/monteCarlo/Makefile.in \
	examples/monteCarlo/src/Makefile.in \
	examples/pyramid/Makefile.in \
	examples/pyramid/src/Makefile.in \
	examples/pyramid/src/physics_1/Makefile.in \
	examples/pyramid/src/physics_1_2/Makefile.in \
	examples/pyramid/src/physics_2/Makefile.in \
	examples/tga/Makefile.in \
	examples/tga/src/Makefile.in

MOSTLYCLEANFILES = $(DX_CLEANFILES)

# Documentation hooks in with the 'html' target
html-local: doxygen-run
doc:        doxygen-run
docs:       doxygen-run

# Required for AX_PROG_DOXYGEN
#include doxygen.am

# Additional files to be deleted by 'make distclean'
DISTCLEANFILES = 


# Required for AX_AM_MACROS
#@INC_AMINCLUDE@
