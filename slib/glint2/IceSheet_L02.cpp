#if 0
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>



typedef boost::numeric::ublas::coordinate_matrix<
	double, boost::numeric::ublas::row_major, 0,
	boost::numeric::ublas::unbounded_array<IceIndex>,
	boost::numeric::ublas::unbounded_array<double>> CooMatrix;

typedef boost::numeric::ublas::coordinate_vector<
	double, 0,
	boost::numeric::ublas::unbounded_array<IceIndex>,
	boost::numeric::ublas::unbounded_array<double>> CooVector;
#endif




// =================================================================


/** Constructs conservative regrid matrix, such that:
     fG = diag(1/weightG) M fI      for each ice sheet
     fI = diag(1/weightI) M^T fG    for each ice sheet

NOTE: This only ADDS to the overall M, scaleI and scaleG, which are
constructed for ALL ice sheets combined.  Indices for THIS ice sheet
start at ix_base.
*/
void IceSheet::ice_xx_exch(
VectorSparseMatrix &M,			// M(X,I)
giss::MapSparseVector<int,double> &weightI,		// 1/x of Ice grid scaling
giss::MapSparseVector<int,double> &weightX)		// 1/x of Exchange grid scaling
const {
	// long const iI_base = i2_to_iI(this->sheetno, 0);
	// long const iX_base = iI_base;
	long const iI_base = 0;
	long const iX_base = 0;

	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		// cell->i = index in atmosphere grid
		long iI = iI_base + cell->j;		// index in ice grid
		long iX = iX_base + cell->index; 	// index in exchange grid
		double val = 1.0;

		M.append_element(iX, iI, val);
		weightI.append_element(iI, val);
		weightI.append_element(iX, val);
	}

}



/** Builds an interpolation matrix to go from height points to ice/exchange grid.
@param dest Controls matrix output to ice or exchange grid

	fG = diag(1/weightG) M fE			# For each ice sheet
	fE = diag(1/weightE) M^T fG			# For each ice sheet

*/
void IceSheet_L0::projelev_xx_iceexch(IceExch dest,
VectorSparseMatrix &M,			// M(G,E)
giss::MapSparseVector<int,double> &weightE,		// Elevation grid scaling
giss::MapSparseVector<int,double> &weightG,		// Ice/Exchange grid scaling
IceExch ice_exch)		// Tells whether G is an ice or exchange grid
{
printf("BEGIN hp_to_iceexch(%s) mask1=%p\n", ice_exch.str(), &*gcm->mask1);
//	long const iE_base = i2_to_iI(this->sheetno, 0);
//	long const iG_base = iE_base;
	long iE_base = 0;		// Never more than one ice sheet at a time in a matrix.
	long iG_base = 0;

	// Offset added to elevation points.  This is to account for the "exta"
	// elevation point that will be provided, if fill_masked.
	// In theory, this parameter can be set independently.
	int ihp_base = (fill_masked ? 1 : 0);

	// ---------------------------------------
	// Handle Z_INTERP or ELEV_CLASS_INTERP


	int nG = niceexch(ice_exch);

	// Interpolate in the vertical
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		int const iA = cell->i;		// GCM (Atmosphere) grid
		int const iI_l = cell->j;		// Ice grid (local to this ice sheet)
		int const iX_l = cell->index;	// Xchange grid (local to this ice sheet)
		int const iG = iG_base + (ice_exch == IceExch::ICE ? iI_l : iX_l);

		if (masked(cell)) {
			if (!fill_masked) continue;

			// This cell is masked out: use special height point 0 (we have ihp_base == 1)
			int ihp = 0;
			int iE = iE_base + gcm->hc_index->ik_to_index(iA, ihp);		// Index in elevation grid

			M(iG,iE) = cell->area;
			weightE(iE) += cell->area;
			weightG(iG) += cell->area;
		} else {
			// This cell not masked: look up elevation point as usual
			double elevation = std::max(elev2(iI_l), 0.0);

			// Interpolate in elevation points
			switch(interp_style.index()) {
				case InterpStyle::Z_INTERP : {
					int ihps[2];
					double whps[2];

					// Get two interpolation weights, and add them to the matrix.
					linterp_1d(gcm->hpdefs, elevation, ihps, whps);
					for (int k=0; k<2; ++k) {
						int iE = iE_base + gcm->hc_index->ik_to_index(i1, ihp_base+ihps[k])
						M(iG,iE) = whps[k] * cell->area;
						weightE(iE) += cell->area;
						weightG(iG) += cell->area;
					}
				} break;
				case InterpStyle::ELEV_CLASS_INTERP : {
					// Add jsut one interpolation point to the matrix.
					int ihps0 = nearest_1d(gcm->hpdefs, elevation);
					int iE = iE_base + gcm->hc_index->ik_to_index(i1, ihps0+hp_base);
					M(iG, iE) = cell->area;
					weightE(iE) += cell->area;
					weightG(iG) += cell->area;
				} break;
			}
		}

	}
}

/**
	fAp = diag(1/weightAp) M fG				# For each ice sheet
	fG = diag(1/weightG) M^t fAp			# For each ice sheet
*/
void IceSheet_L0::iceexch_xx_projatm(
VectorSparseMatrix &M,				// M(A,G)
giss::MapSparseVector<int,double> &weightG,		// Ice/Exchange grid scaling
giss::MapSparseVector<int,double> &weightAp,		// Projected atmosphere grid scaling
IceExch ice_exch)		// Tells whether G is an ice or exchange grid
{
printf("BEGIN IceSheet_L0::ice_to_projatm %ld %ld\n", n1(), n4());

	// long const iG_base = i2_to_iI(this->sheetno, 0);
	// long const iAp_base = iG_base;
	long const iG_base = 0;
	long const iAp_base = 0;

	// ============= ice_to_atm (with area1 scaling factor)
	// Area-weighted remapping from exchange to atmosphere grid is equal
	// to scaled version of overlap matrix.

	int nG = niceexch(ice_exch);

	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		// Exchange Grid is in Cartesian coordinates
		// cell->i = index in atmosphere grid
		// cell->j = index in ice grid
		// cell->index = index in exchange grid
		const int iG = iG_base + (ice_exch == IceExch::ICE ? cell->j : cell->index);
		const int iAp = iAp_base + cell->i;

		M.add(iAp,iG, cell->area);
		weightG(iG) += cell->area;
		weightAp(iAp) += cell->area;
	}
}

