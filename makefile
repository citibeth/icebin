# This branch of HTML documentation should be checked out in relation to
# the main branch:
#    dir/glint2 = main branch
#    dir/glint2-html = gh-pages branch (see GitHub for gh-pages info)
#
# See:
#    http://rickfoosusa.blogspot.com/2011/10/howto-use-doxygen-with-github.html
#    http://docmartinscodingstuff.blogspot.com/2012/09/github-pages-doxygen-minor-snag.html

doxygen-doc :
	\rm -rf ../glint2/doxygen-doc
	cd ../glint2; make doxygen-run
	\rm -rf doxygen
	mv ../glint2/doxygen-doc/html doxygen
	mv ../glint2/doxygen-doc/glint2.tag .

add-doxygen :
	cd doxygen; git add -u
	git add `find doxygen -name '*'`

push :
	git push origin gh-pages
