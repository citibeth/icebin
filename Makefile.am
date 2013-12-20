# this is example-file: Makefile.am

# the subdirectories of the project to go into
SUBDIRS = 	slib \
			sbin \
			pyext

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = glint2.pc

include am_doxygen.am

if DX_COND_man

# You'd probably want to post-process man pages and installed the patched
# versions.
dist_man1_MANS = @DX_DOCDIR@/man/man1/doxample.c.1

$(dist_man1_MANS): doxygen-doc

endif

MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = cleanup bootstrap $(TESTS) $(DX_CONFIG)

modele_const.py : parse_fconst.py
	python parse_fconst.py >modele_const.py
