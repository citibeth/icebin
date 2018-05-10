.. _sparse_matrices

Sparse Sparse Matrices
======================

The matrices used and produced by *IceBin* are sparse in two ways:
they are sparse matrices in the traditional sense, and they also use
sparse indexing.

Sparse Matrices
---------------

A *sparse matrix* is one in which most of its elements are zero.
Sparse matrices are occur commonly in physics problems; for example,
:math:`\matrix{O}^{BA}` is sparse.  More precisely, a matrix of rank
:math:`n \times m` is sparse if asymptotically, the number of non-zero
elements is less than :math:`\Theta(mn)`.  This means that not only
are most of the elements zero, but the *proportion* of non-zero
elements goes down as larger versions of the matrix are realized.  For
problems that produce sparse matrices, it is therefore *essential* to
use proper sparse matrix representations when solving them.

Sparse Indexing
---------------

Not all indices ("grid cells") of a vector space are always used.  For
example, a grid might be created to represent ice cover globally.
Since only 12% of the Earth is covered in ice, only 12% of grid cells
in this grid will be used.  Similarly, suppose that 20 elevation
classes are defined for ice-covered GCM grid cells.  Typically at most
2-5 elevation (out of 20) classes might be turned on for any
particular GCM grid cell.  Taking both these effects into account,
only about 1% of elevation classes are used in a typical ModelE run;
the rest have an index and may have memory allocated to them, but are
not involved in any regridding computation.  More abstractly, they are
part of every regridding matrice's nullspace.

Unlike sparse matrices, the sparseness of sparse indexing does not
change with scale.  If 12% of the earth is ice-covered, then
approximately 12% of grid cells will be used, no matter the resolution
of those grid cells.

IceBin
------

*IceBin* is made to construct sparse matrices on sparsely indexed
vector spaces, such as ice or elevation grids.  Each vector space
("grid") is mapped internally to a dense subspace containing only the
dimensions ("grid cells") that are in use.  *IceBin* maintains a
mapping between the original ("sparse") and dense version of each
vector space.  Depending on the situation, matrices are converted back
to sparse indexing when returned to the user; or they are retunred
"as-is", along with the sparse-to-dense mapping.


ASCII Notation
--------------

* The varaiable ``dimA`` maps between sparse and dense indexing for vector space :math:`\cal{A}`.
* The suffix ``_s`` labels a variable as using sparse (native) indexing.
* The suffix ``_d`` labels a varaible as using dense (mapped) indexing.

