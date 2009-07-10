## Process this file with automake to produce Makefile.in

# Don't require GNU-style project
AUTOMAKE_OPTIONS = foreign

# Build in these directories:
SUBDIRS = libs examples test/t01_valid_cycle/valid_cycle test

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/missing

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES =     \
	aclocal.m4             \
	aminclude.am           \
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

# Additional files to be deleted by 'make distclean'
DISTCLEANFILES = \
	aminclude.am