.. _matrix_formats

Matrix Formats
==============

*IceBin* generates matrices and vectors that are store on-disk in one of
two formats: *eigen* or *compressed*.

Eigen Format
------------

Eigen Format stores matrices in NetCDF in a straightforward,
flexible format.  Consider the following NetCDF file, for
example, which defines a matrix ``BvA.M`` and vectors ``BvA.Mw`` and
``BvA.wM``.

.. code-block:: none

   netcdf __linear1 {
   dimensions:
           dimB.dense_extent = 5 ;
           dimA.dense_extent = 6 ;
           BvA.M.nnz = 6 ;
           BvA.M.rank = 2 ;
   variables:
           int64 dimB(dimB.dense_extent) ;
                   dimB:sparse_extent = 40LL ;
           int64 dimA(dimA.dense_extent) ;
                   dimA:sparse_extent = 50LL ;
           int BvA.info ;
                   BvA.info:type = "EIGEN" ;
                   BvA.info:conservative = 0 ;
                   string BvA.info:dim_names = "BvA.dimB", "BvA.dimA" ;
           int64 BvA.M.info ;
                   BvA.M.info:shape = 5LL, 6LL ;
                   BvA.M.info:conservative = "f" ;
           int BvA.M.indices(BvA.M.nnz, BvA.M.rank) ;
           double BvA.M.values(BvA.M.nnz) ;
           double BvA.Mw(BvA.dimA.dense_extent) ;
           double BvA.wM(BvA.dimB.dense_extent) ;
   }

In this case, two vector spaces are involved: :math:`\cal{A}` and
:math:`\cal{B}`.  The variable ``dimA`` defines a "dense" vector space
by listing the indices of :math:`\cal{A}` that are used in vectors and
matrices.  The attribute ``dimA:sparse_extent`` give the rank of the
original sparse vector space :math:`\cal{A}`; and the dimension
``dimA.dense_extent`` gives the rank of the densified vector space
defined by ``dimA``.

Vectors are stored in regular dense format in the dense space
``dimA``.  For example, if ``BvA.Mw[0] == 28.9`` and ``dimA[0] ==
17``, then the element at index ``17`` in the original sparse vector
space :math:`\cal{A}` equals ``28.9``.

Matrices (tensors) are stored in as coordinate-format sparse matrices
within the densified vector space ``dimA`` and ``dimB``.  For the
matrix ``BvA.M``, coordinate indices are stored in ``BvA.M.indices``
and the sparse matrix value at those indices is stored at
``BvA.M.indices``.  In addition, the attribute
``BvA.M.info:conservative`` declares whether this (regridding) matrix
is conservative; non-conservative regridding matrices require a
conservation correction when applied.

Compressed Format
-----------------

Compressed format uses `zlib <https://www.zlib.net/>`_ compression to
reduce the on-disk and in-memory footprint of large matrices, while
still allowing for them to be multiplied by vectors.  Compressed
format vectors and matrices are stored in the original sparse vector
space :math:`cal{A}`, allowing for easy use in applications where
sparse vector spaces are not densified (such as ModelE).  Here is an example:

.. code-block:: none

   netcdf __linear3 {
   dimensions:
           BvA.wM.indices.zsize = 29 ;
           BvA.wM.values.zsize = 36 ;
           BvA.M.indices.zsize = 41 ;
           BvA.M.values.zsize = 30 ;
           BvA.Mw.indices.zsize = 30 ;
           BvA.Mw.values.zsize = 33 ;
   variables:
           int BvA.info ;
                   BvA.info:type = "COMPRESSED" ;
                   BvA.info:conservative = 0 ;
           int BvA.wM.info ;
                   BvA.wM.info:format = "ZARRAY" ;
                   BvA.wM.info:rank = 1 ;
                   BvA.wM.info:nnz = 3LL ;
                   BvA.wM.info:shape = 40LL ;
           ubyte BvA.wM.indices(BvA.wM.indices.zsize) ;
           ubyte BvA.wM.values(BvA.wM.values.zsize) ;
           int BvA.M.info ;
                   BvA.M.info:format = "ZARRAY" ;
                   BvA.M.info:rank = 2 ;
                   BvA.M.info:nnz = 6LL ;
                   BvA.M.info:shape = 40LL, 50LL ;
           ubyte BvA.M.indices(BvA.M.indices.zsize) ;
           ubyte BvA.M.values(BvA.M.values.zsize) ;
           int BvA.Mw.info ;
                   BvA.Mw.info:format = "ZARRAY" ;
                   BvA.Mw.info:rank = 1 ;
                   BvA.Mw.info:nnz = 4LL ;
                   BvA.Mw.info:shape = 50LL ;
           ubyte BvA.Mw.indices(BvA.Mw.indices.zsize) ;
           ubyte BvA.Mw.values(BvA.Mw.values.zsize) ;
   }

Vectors and matrices are both stored as *compressed tensors* of rank 1
or 2, respectively.  In the above case, the variables
``BvA.M.indices`` and ``BvA.M.values`` store the matrix in coordinate
format; however, the contents of those arrays is not directly readable
without first running through decompression.  *IceBin* provides
libraries, accessible from Python and C++, that will decompress these
sparse matrices as needed.


Python API
----------

The following sample program demonstrates how to read and use sparse
matrices to NetCDF using Python.  The Python API for both of matrices
is the same, shielding the user from implementation differences
between the two.

*IceBin* typically generates a regridding matrix, along with *weight
vectors* for the two dimensions involved in the matrix.  Depending on
the case, the weight vectors may or may not be the same as the sum of
rows and columns of the matrix.  The Python class
``ibmisc.linear_Weighted`` holds a matrix plus two vectors, and allows
basic multiplication and dot product operations on them.

Loading
^^^^^^^

To load a ``linear_Weighted`` from a NetCDF file, use (for example):

.. code-block:: python

   with ibmisc.NcIO('file.nc', 'r') as ncio:
       lw = ibmisc.nc_read_weighted(ncio, 'BvA')

Dot Products
^^^^^^^^^^^^

The matrix in ``linear_Weighted`` may be applied to a vector using the
``apply_M()`` method.  Note that the input ``A_s`` is in the
*original* sparse vector space :math:`\cal{A}`:


.. code-block:: python

   def apply_M(self, A_s, fill=np.nan, bool force_conservation=True):
       """Applies the regrid matrix to A_s.
       A_s: Either:
           - A single vector (1-D array) to be transformed.
           - A 2-D array of row vectors to be transformed.
       fill:
           Un-set indices in output array will get this value.
       force_conservation: bool
           If M is not conservative, apply a conservation correction
           at the end"""


Similarly, the method ``apply_weight()`` takes inner products with the
weight vectors:

.. code-block:: python

   def apply_weight(self, int dim, A_s):

       """Computes dot product of a weight vector with A_s.
       dim:
           0 = use weight vector for B (output) dimension
           1 = use weight vector for A (input) dimension
       A_s: Either:
           - A single vector (1-D array) to be transformed.
           - A 2-D array of row vectors to be transformed."""

   def apply_wM(self, A_s):
       return self.apply_weight(0,A_s)
   def apply_Mw(self, A_s):
       return self.apply_weight(1,A_s)


Extraction
^^^^^^^^^^

The method ``to_coo()`` returns the matrix in standard Python
coordinate format, as a ``scipy.sparse.coo_matrix`` object.
Similarly, the method ``get_weights(dim)`` returns the weight vector
(0=B, 1=A) as a ``numpy.array``.

Example
^^^^^^^

See `here
<https://github.com/citibeth/ibmisc/blob/master/pylib/tests/test_cython.py>`_
for a working example of Python code that uses this API.
