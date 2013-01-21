# nobase_include_HEADERS = mylib/myClass.hpp
# nobase_include_HEADERS = `find . -name '*.hpp'`

nobase_include_HEADERS = $(shell \
		find -E . -regex '.*\.(hpp|h|mod)' \
	)
