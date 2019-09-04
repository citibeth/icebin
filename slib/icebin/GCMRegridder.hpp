/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <functional>
#include <unordered_set>
#include <blitz/array.h>

#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/IndexSet.hpp>
#include <spsparse/eigen.hpp>
#include <ibmisc/linear/tuple.hpp>

#include <icebin/IceRegridder.hpp>
#include <icebin/RegridMatrices.hpp>
#include <icebin/RegridMatrices_Dynamic.hpp>

/** 

This module produces regridding matrices between the following vector
spaces (grids):

     Description                    Geometry
---------------------------------------------
A  = Atmosphere grid                sphere
Ap = Projected atmosphere grid      plane
E  = Elevation grid                 sphere
Ep = Projected elevation grid       plane
I  = Ice grid                       plane

The atmosphere grid ``A`` is assumed to be L0 (aribtrary polygonal
grid cells with constant value inside).  However, ``I`` may use any
set of basis functions; for example, a finite element mesh.  The
following "grids" are also used sometimes:

G = Interpolation Grid
    --> For constant-value (L0) ice grids, this is either the ice grid itself,
    or the exchange grid (overlap) between ``I`` and ``A``.


Regridding Basics
-----------------

* The symbol ``nX = |X|`` is the number of basis functions (grid cells) in
  a given vector space.  Not all of these need to be *realized*: an
  ice grid, for example, might have some cells that are masked out
  because they don't contain ice.

* A regridding matrix named ``BvA`` (size ``nb`` * ``nA``) regrides
  vectors from grid ``A`` to ``B``.

* ``wA`` is the vector (or diagonal matrix) containing the weight
  (area) of each basis function (grid cell) of ``A``.

* ``sA`` is the inverse of ``wA``... it is a "scaling" matrix, not a
  "weight" matrix.

* All regridding from ``A`` to ``B`` is accomplished by taking the
  inner product of the ``B`` basis functions with the `f(x,y)` defined
  by ``A``.  This is accomplished by computing the unscaled matrix
  ``BvA``, where the ``(j,i)`` element is the inner product of the
  ``j``th basis function of ``B`` with the ``i``th basis function of
  ``A``.  This has some interesting properties:

  #. The unscaled ``AvB`` and unscaled ``BvA`` are transposes of each
     other.  Computing the unscaled ``AvB`` there fore the first step
     to computing the final scaled ``AvB`` or ``BvA``.

  #. This matrix produces *unscaled* results; that is, if the vector
    ``a`` is in [kg m-2], then the product ``BvA * a`` will be in
    [kg].  One must divide by ``wB`` to get a scaled result.

  #. The scaling weights for an unscaled ``BvA`` may be obtained by
     summing over the rows or columns.  This is almost always
     preferable to using overall weights ``wA`` and ``wB`` of the grid
     cells because it properly accounts for cells of ``A`` and ``B``
     that do not overlap the other grid.  We call these (diagonal)
     weight matrices ``wAvB`` and ``wBvA``.  Their inverses are
     ``sAvB`` and ``sBvA``.

* Putting it together... if one has an unscaled ``BvA`` and have
  summed its rows to obtain ``wBvA``, then the scaled regridding
  matrix will be ``sBvA * BvA``.

* The atmosphere and elevation grids ``A`` and ``E`` exist on a
  sphere, not a plane.  They therefore cannot interact with the ice
  grid, which exists on a plane.  This issue is resolved by
  *projecting* ``A`` and ``E`` onto the plane.  The result is the
  grids ``Ap`` and ``Ep``, which we *can* regrid to/from ``I``.  A
  diagonal scaling matrix ``sApvA`` converts from ``A`` to ``Ap`` (or
  ``E`` to ``Ep``).  There is symmetry here: note that ``sApvA`` =
  ``wAvAp``

Regrid Hierarchy
----------------

Regridding matrices are computed in a three-level hierarhcy:

* ``IceRegridder``: Computes a variety of basic unscaled regrid
  component matrices: ``sApvA``, ``GvAp``, ``sEpvE``, ``GvEp``,
  ``GvI``.  These are computed based on detailed knowledge of the
  grids; they will vary depending on the grid's parameterization,
  outlines, etc.  Note the symmetry between ``A`` and ``E`` here:
  although the details of how matrices regridding to ``A`` and ``E``
  may vary, they are *used* in algebraically similar ways.

* Ur ("original") Matrices: The ``UrAE`` class encapsulates the algebraic symmetry
  between ``A`` and ``E``.  Depending on how it's instantiated, it
  will provide *either* ``GvAp`` and ``sApvA`` *or* ``GvEp`` and
  ``sEpvE`` (taken from the lower-level functions in the step above).

  The Ur instance is passed into higher-level functions that will
  produce regridding matrices going to/from either ``A`` or ``E``.
  For simplicity sake, the grid is always called ``A`` in those
  higher-level routines.

* Matrix Combiner Functions: Three free functions are responsibile for
  multiplying togehter lower-level matrices to compute final regrid *
  matrices.  The final matrices may be either scaled or unscaled; and
  have project issues (``Ap`` vs ``A``) corrected or not corrected.
  Due to the symmetry between ``A`` and ``E``, each one does double
  duty:

      Function               Final Matrices Computed
      ------------------------------------------------
      compute_AEvI()         AvI, EvI
      compute_IvAE()         IvA, IvE
      compute_EvA()          EvA, AvE

* ``RegridMatrices``: The ``RegridMatrices`` class wraps all this
  structure in a simple, top-level interface to generate regridding
  matrices, ``RegridMatrix::regrid()``.  For example:

  .. code-block:: c++

     GCMRegridder gcm_regridder(...);
     RegridMatrices rm(gcm_regridder.sheet("greenland"));
     // rm.regrid("AvI", scale=true, correctA=true)
     auto AvI(rm.regrid("AvI", true, true));
     AvI.M        // The matrix
     AvI.weight   // The weight vector


Sparse Computations
-------------------

Matrices and vectors are represented in coordinate fashion by
``spsparse::SparseMatrix`` and ``spsparse::SparseVector``.  This
representation is good for two reasons:

  1. It is a convenient and flexible way to manage sparse matrices.

  2. It allows matrices to be sparse in their index set, as well as in
     rows and columns: the size of each dimension of the matrix may be
     arbitrarily large, as long as the number of non-zero elements
     fits in memory.  This is useful for implementing aribtrary grids,
     where grid cell indices may become large.

When multiplying these sparse indices and vectors, they are first
converted into Eigen sparse matrices --- which require a dense index
set so that a row or column may be represented as a dense array.
Native indices are mapped to a dense set for Eigen computations, and
then mapped back when the answer is converted back to a ``spsparse``
structure.

*/


namespace icebin {

class GCMRegridder;

// ----------------------------------------------------

extern ibmisc::Indexing derive_indexingE(
    ibmisc::Indexing const &indexingA,      // im,jm,...
    ibmisc::Indexing const &indexingHC);    // iA,iHC


// ----------------------------------------------------

/** Used to index arrays that are done for A and E grids */
namespace GridAE {
    static const int A = 0;
    static const int E = 1;
    static const int count = 2;
}

/** Generates the matrices required in the GCM */
class GCMRegridder
{

public:
    /** Only used when first creating a GCMRegridder (Cython).
    This field is not stored or loaded with ncio(). */
//    std::unique_ptr<Grid> fgridA;

    AbbrGrid agridA;

//  ibmisc::Domain domainA;                // What's in our MPI halo?

    /** Should we correct for projection and geometric error?
    NOTE: This is read out of the IceBin file on disk, but not used
          when generating regrid matrices (or anywhere else).  It may be
          queried by an application program wishing to know the value
          of correctA stored in the file; and then sent into the regridding
          routines. */
    bool correctA;

    /** Each grid cell iE in the elevation grid may be thought of as
        related to a grid cell iA in the atmosphere grid, along with
        an elevation class iHP.  This converts converts indices
        between (iA, iHP) <--> (iE).  The iA index may be broken down
        further into grid-specific coordinates using gridA->indexing.

        (iA, iHP) zero-based
    */
    ibmisc::Indexing indexingHC;

    // Derived from gridA.indexing and indexingHC
    // (im,jm,...,ihp) zero-based indices
    ibmisc::Indexing indexingE;

    /** Selection function so aid when A/E stuff uses the same code,
        indexed by GridAE::A or ::E. */
    ibmisc::Indexing const &indexing(int iAE) const
        { return (iAE == GridAE::A ? agridA.indexing : indexingE); }

protected:
    /** Ice sheets stored by index defined in sheets_index */
    ibmisc::IndexedVector<std::string, std::unique_ptr<IceRegridder>> *_ice_regridders;

public:

    /** Position of height points in elevation space (same for all GCM
    grid cells and all ice sheets) */
    std::vector<double> _hcdefs; // [nhc]

    std::vector<double> const &hcdefs() const
        { return _hcdefs; }

    ibmisc::IndexedVector<std::string, std::unique_ptr<IceRegridder>> &ice_regridders()
        { return *_ice_regridders; }
    ibmisc::IndexedVector<std::string, std::unique_ptr<IceRegridder>> const &ice_regridders() const
        { return *_ice_regridders; }

    virtual ~GCMRegridder();

    /** @return Number of elevation points for grid cells in general */
    /** @return Number of elevation points for a given grid cell */
    unsigned int nhc(int i1) const { return (unsigned int)_hcdefs.size(); }
    unsigned int nhc() const { return nhc(-1); }

    unsigned long nA() const { return agridA.dim.sparse_extent(); }
    unsigned long nE() const { return nA() * nhc(-1); }
    size_t nI(int sheet_index) const { return ice_regridders()[sheet_index]->nI(); }

    template<class AccumT>
    void wA(AccumT &&accum, std::string const &ice_sheet_name, bool native);

    // ==================== Compute merged matrices for all ice sheets
    // "Extra" operations; such as adding legacy ice, ECs (for legacy ice), etc.
    // can happen here.


    /** Produce regridding matrices for this setup.
    Do not change elevmaskI to dense indexing.  That would require a
    SparseSet dim variable, plus a dense-indexed elevation.  In the end,
    too much complication and might not even save RAM. */
    virtual std::unique_ptr<RegridMatrices_Dynamic> regrid_matrices(
        int sheet_index,
        blitz::Array<double,1> const &elevmaskI,
        RegridParams const &params = RegridParams()) const = 0;

    /**
    @param rw_full If true, read the entire data structure.  If false (i.e. we
                   are using MPI and this is not the root), then avoid reading
                   grid details, etc.
    */
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};

template<class AccumT>
void GCMRegridder::wA(AccumT &&accum, std::string const &ice_sheet_name, bool native)
{
    IceRegridder *ice = &*ice_regridders().at(ice_sheet_name);

    auto &areas(native ? agridA.native_area : ice->gridA_proj_area);
    for (int id=0; id<agridA.dim.dense_extent(); ++id) {
        auto index = agridA.dim.to_sparse(id);
        accum.add({index}, areas(id));
    }
}

// ----------------------------------------------------------------------
/** Generates the matrices required in the GCM */
class GCMRegridder_Standard : public GCMRegridder
{

    /** Creates an (index, name) correspondence for ice sheets. */
    ibmisc::IndexSet<std::string> mem_sheets_index;

    /** Ice sheets stored by index defined in sheets_index */
    ibmisc::IndexedVector<std::string, std::unique_ptr<IceRegridder>> mem_ice_regridders;

public:

    /** Constructs a blank GCMRegridder.  Typically one will use
        ncio() afterwards to read from a file. */
    GCMRegridder_Standard() :
        mem_ice_regridders(mem_sheets_index)
    {
        _ice_regridders = &mem_ice_regridders;
    }

    void clear();

    /** Used to initially construct a GCMRegridder.  Most applications
        programs will not use this; they will use ncio() to read one
        in from a file.

    @param _gridA
        The atmosphere grid to use
    @param _hcdefs
        Elevation of each elevation class to use.  (The
        same definitions will be used for every grid cell).
    @param _indexingHC
        How to convert (iA, iHP) <--> (iE) */
    void init(
        AbbrGrid &&_agridA,
//        std::unique_ptr<Grid> &&_gridA,
//      ibmisc::Domain &&_domainA,     // Tells us which cells in gridA to keep...
        std::vector<double> &&_hcdefs,
        ibmisc::Indexing &&_indexingHC,
        bool _correctA);

    // -----------------------------------------
    // These subroutines used for initial construction...

    void add_sheet(std::unique_ptr<IceRegridder> &&regridder)
    {
        regridder->gcm = this;
        size_t ix = ice_regridders().index.insert(regridder->name());
        ice_regridders().push_back(std::move(regridder));
    }

    /** Convenience method, for when the IceRegridder does not yet
        have a name. */
    void add_sheet(std::string name, std::unique_ptr<IceRegridder> &&regridder)
    {
        regridder->_name = name;
        add_sheet(std::move(regridder));
    }

    /** Produce regridding matrices for this setup.
    Do not change elevmaskI to dense indexing.  That would require a
    SparseSet dim variable, plus a dense-indexed elevation.  In the end,
    too much complication and might not even save RAM. */
    std::unique_ptr<RegridMatrices_Dynamic> regrid_matrices(
        int sheet_index,
        blitz::Array<double,1> const &elevmaskI,
        RegridParams const &params) const;

    /** Removes unnecessary cells from the A grid
    @param keepA(iA):
        Function returns true for cells we wish to keep. */
    void filter_cellsA(std::function<bool(long)> const &keepA);

    /** Higher-level filter, filter out cells not in our MPI domain.
    @param domainA Description of our MPI domain.
    Indices are in C order with 0-based indexing. */
    void filter_cellsA(ibmisc::Domain const &domainA);

    typedef ibmisc::DerefRandomAccessIter<
        const IceRegridder,
        typename std::vector<std::unique_ptr<IceRegridder>>::const_iterator
    > const_iterator;

    // -----------------------------------------
    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};  // class GCMRegridder_Standard
// ===========================================================



}   // namespace
