## Process this file with automake to produce Makefile.in

# Don't require GNU-style project
AUTOMAKE_OPTIONS = foreign

# Build in these directories:
SUBDIRS = libs examples

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
	libs/queso/Makefile.in \
	libs/queso/src/Makefile.in \
	libs/queso/inc/Makefile.in \
	examples/Makefile.in   \
	examples/queso/Makefile.in \
	examples/queso/normal/Makefile.in \
	examples/queso/normal/src/Makefile.in \
	examples/queso/pyramid/Makefile.in \
	examples/queso/pyramid/src/Makefile.in \
	examples/queso/pyramid/src/physics_1/Makefile.in \
	examples/queso/pyramid/src/physics_1_2/Makefile.in \
	examples/queso/pyramid/src/physics_2/Makefile.in \
	examples/queso/tga/Makefile.in \
	examples/queso/tga/src/Makefile.in

# Additional files to be deleted by 'make distclean'
DISTCLEANFILES = \
	aminclude.am