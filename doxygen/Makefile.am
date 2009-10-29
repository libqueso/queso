## Process this file with automake to produce Makefile.in

# Don't require GNU-style project
AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = queso.page txt_common fig_common

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

