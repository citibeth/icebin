#!/bin/env python3
#
import sys
import os

# http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
def search_file(filename, search_path):
     """Given a search path, find file
     """
     file_found = 0
     paths = search_path.split(os.pathsep)
     for path in paths:
         if os.path.exists(os.path.join(path, filename)):
             file_found = 1
             break
     if file_found:
         return os.path.abspath(os.path.join(path, filename))
     else:
         return None


fname = search_file(sys.argv[1], os.environ[sys.argv[2]])
if fname is None:
    sys.exit(-1)
else:
    print(fname)

