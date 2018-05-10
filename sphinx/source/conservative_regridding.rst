.. _conservative_regridding


Generalized Conservative Regridding
===================================

Regridding problem is the problem of transforming a vector
:math:`\vec{f^A}` in vector space :math:`\cal{A}` into an "equivalent"
vector :math:`\vec{f^B}` in vector space :math:`\cal{B}`.  Because
:math:`\cal{A}` is a vector spaces, it has dimensionality
:math:`|\cal{A}|` and basis functions, denoted as
:math:`a_i(\vec{x})` for :math:`1 \leq x \leq |\cal{A}|`.

.. note::
   * Regridding transformations are almost always *linear*, and therefore
     representable by a matrix.  All regridding operations we consider
     here will be linear.

   * Practitioners typically think about regridding in terms of
     *grids*.  Grids can be more formally represented by the more
     general notion of a *vector space*, where each "grid cell" is
     represented by a basis function.

The notion of "equivalence" may be formalized by insisting that
:math:`\cal{N}^A(\vec{f^A}) = \cal{N}^B(\vec{f^B})` where
:math:`\cal{N}^A` and :math:`\cal{N}^B` are norms on the vector
spaces.  We can define the vector of weigths for the basis functions:

.. math::
   \vec{|a|} = |a_i| = \int_A a_i(\vec{x}) dA

Using this, we can define a norm that corresponds to our idea of
"conservative", i.e. it measures the amount of "stuff" represented by
a vector (assuming :math:`\vec{f^A}` uses intrusive units).

.. math::
   \cal{N}^A(\vec{f^A}) = \vec{|a|} \cdot \vec{f^A}

.. note::

   Intrusive units are units, such as *kg/m^2*, express the amount of
   a substance in a grid cell by the *quantity* of substance, divided
   by the area/volume.  Climate models typically use intrusive units
   for all quantities.

Overlap Matrix
--------------

We define an *overlap matrix* :math:`\matrix{O}^{BA}`, or *unscaled
regridding matrix* to be:

.. math::
   \matrix{O}^{BA}_{ji} = b_j(\vec{x}) \cdot a_i(\vec{x}) = \int_A b_j(\vec{x}) a_i(\vec{x}) dA

We define *weight vectors* for vector spaces :math:`\cal{A}` and :math:`\cal{B}` as:

.. math::
   \vec{w^B} = \matrix{O}^{BA} \cdot \vec{1} \\
   \vec{w^A} = \vec{1} \cdot \matrix{O}^{BA}

:math:`\vec{w^B}` can be interpreted as the area of grid cells
in :math:`\cal{B}` that overlap :math:`\cal{A}`.  If the entire grid
:math:`\cal{B}` overlaps :math:`\cal{A}`, then :math:`\vec{w^B} =
\vec{b}`.  However, this is often not the case: for example, an ice
sheet may be define on only a portion of a grid, with the other grid
cells unused.  The effective grid consists of only ice sheet-covered
cells, and :math:`\vec{w^B}` must be used instead of :math:`\vec{b}`
in evaluating conservation properties on this grid.

Scaled Regridding Matrix
------------------------

We define *scale vectors* as the inverse of weight vectors
:math:`s^A_i = 1 / w^A_i`.  A conservative (or "scaled") regridding
matrix may be computed as:

.. math::
   \matrix{M}^{BA} = \textit{diag}(s^B) \matrix{O}^{BA}

It can be used to regrid:

.. math::
   \vec{f^B} = \matrix{M}^{BA} \vec{f^A}


with the property:

.. math::
   \vec{b} \cdot \vec{f^B} = \vec{a} \cdot \vec{f^A}

It is important to note that *all* conservative regridding matrices
published in the literature follow this pattern.  The only
"interesting part" is the computation of the inner products
:math:`\int_A b_j(\vec{x}) a_i(\vec{x}) dA`.

Regridding Areas vs. Quantities
-------------------------------

The above formulas work if :math:`\vec{f^A}` uses intrusive units; eg
*kg/m^2*; and they can also be trivially adapted to the case where
:math:`\vec{f^A}` is not scaled by area, eg. in units of *kg*.  If :math:`\vec{f^A}` is in units of area (*m^2*), then the following formula applies instead:

.. math::
   \vec{f^B} = \matrix{O}^{BA} \textit{diag}(\vec{s^A}) \vec{f^A}

Symmetries and Identities
-------------------------

Regridding matrices have a fundamental symmetry based on the fact
that

.. math::
   \matrix{O}^{AB} = (\matrix{O}^{BA})^T

Thus, :math:`\matrix{M}^{BA}` and :math:`\matrix{M}^{AB}`
may both be computed from :math:`\matrix{O}^{BA}` and its weight vectors
:math:`w^B` and :math:`w^A`.  This suggests a layered system to
compute scaled regridding matrices: core subroutines that compute
:math:`\matrix{O}^{BA}`, and then outer subroutines that scale and (possible)
transpose depending on whether :math:`\matrix{M}^{BA}` or :math:`\matrix{M}^{AB}` is
required.

ASCII Notation
--------------

The following notation is used for variables in *IceBin* code:

* ``BvA`` or ``BvA.M`` = :math:`\matrix{O}^{BA}` or :math:`\matrix{M}^{BA}`, depending on context.
* ``wB`` = :math:`\vec{w^B}`
* ``sB`` = :math:`\vec{s^B}`
* ``wBvA`` or ``BvA.wM`` = :math:`\vec{w^B}` computed from the matrix :math:`\matrix{O}^{BA}`.
* ``sBvA`` = :math:`\vec{s^B}` computed from the matrix :math:`\matrix{O}^{BA}`.
* ``BvAw`` or ``BvA.Mw`` = :math:`\vec{w^A}` computed from the matrix :math:`\matrix{O}^{BA}`.
* ``BvAs`` = :math:`\vec{s^A}` computed from the matrix :math:`\matrix{O}^{BA}`.

In the code, final regridding matrices are typically assembled by
computing lower-level matrices and then multiplying them together.
The multiplication is done with the *Eigen* package, thereby
elucidating the core mathematical nature of what is being multiplied.

Regridding Chains
-----------------

In some cases, a desired regridding matrix is not available directly between two vector spaces.  Consider the following example, with grids :math:`\cal{A}`, :math:`\cal{B}` and :math:`\cal{C}`.  In this case, we know how to compute :math:`\matrix{O}^{BA}` and :math:`\matrix{O}^{CB}`:

.. http://graphs.grevian.org/example
.. http://graphs.grevian.org/reference

.. graphviz::

   digraph {
      rankdir=RL;
      A -> B[label="BvA"]
      B -> C[label="CvB"]
   }

Linearity of these transformations allows us to compute the scaled regrid matrix:

.. math::
   \matrix{M}^{CA} = \mbox{sCvB} \cdot \mbox{CvB} \cdot \mbox{sBvA} \cdot \mbox{BvA}

The equivalent unscaled version can be obtained by omitting the last term on the left:

.. math::
   \matrix{O}^{CA} = \mbox{CvB} \cdot \mbox{sBvA} \cdot \mbox{BvA}

Regridding matrices can therefore be assembled based on graph
diagrams, like the one above, showing the vector spaces and the
available base-level regridding matrices we know how to compute.
