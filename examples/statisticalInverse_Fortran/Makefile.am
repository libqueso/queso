## Process this file with automake to produce Makefile.in

#Build in these directories:
SUBDIRS    = src
EXTRA_DIST = tests

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

