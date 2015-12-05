



/** Adds to a sparse matrix that projects from atmosphere to projected
atmosphere.  Note there is one projected atmosphere per ice sheet. */

void atm_project(			// Transformation: Ap --> A
VectorSparseMatrix &M,					
ProjCorrect direction,
// ---- Only used if direction == ProjCorrect::PROJ_TO_NATIVE
giss::MapSparseVector<int,double> const *iweightAp,		// IN : Weight vector from iceexch_xx_projatm: portion of this projected atm cell used
giss::MapSparseVector<int,double> *weightA)				// OUT: Scaling to use when applying this transformation
{
	// long const iAp_base = i2_to_iI(this->sheetno, 0);
	long const iAp_base = 0;		// Never more than one ice sheet at a time in a matrix.

	giss::Proj2 proj;
	gcm->grid1->get_ll_to_xy(proj, grid2->sproj);

	// TODO: Consider copying iweightAp to a hash_map<int,double> for greater speed in loop.

	for (auto cell = gcm->grid1->cells_begin(); cell != gcm->grid1->cells_end(); ++cell) {
		if (gcm->masked1(cell->index)) continue;

		double native_area = cell->area;
		double proj_area = area_of_proj_polygon(*cell, proj);

		long iA = cell->index;
		long iAp = iAp_base + iA;


		if (direction == ProjCorrect::NATIVE_TO_PROJ) {
			double coeff = native_area / proj_area;
			M.add(iAp, iA, coeff);
			// No weight vector needed here because it would be trivial.
			// In theory, the would match what's below, but things cancel out:
			// double coeff = weightA(iA) * native_area / proj_area;
			// M.add(iAp, iA, coeff);
			// weightA.add(iAp, weightA(iA));    // Only one thing gets added here, so no need for weight vector!

		} else {
			double coeff = iweightAp(iAp) * proj_area / native_area;
			M.add(iA, iAp, coeff);
			// After dividing by weightA, mean of final matrix will be (proj_area / native_area)
			weightA.add(iA, iweightAp(iAp));
		}
	}
}





/** Adds to a sparse matrix that projects from atmosphere to projected
atmosphere.  Note there is one projected atmosphere per ice sheet. */

void elev_project(			// Transformation: Ep --> E
VectorSparseMatrix &M,					
ProjCorrect direction,
// ---- Only used if direction == ProjCorrect::PROJ_TO_NATIVE
giss::MapSparseVector<int,double> const *iweightEp,		// IN : Weight vector from projelev_xx_iceexch: portion of this projected atm cell used
giss::MapSparseVector<int,double> *weightE)				// OUT: Scaling to use when applying this transformation
{
	// long const iEp_base = i2_to_iI(this->sheetno, 0);
	long const iEp_base = 0;		// Never more than one ice sheet at a time in a matrix.

	giss::Proj2 proj;
	gcm->grid1->get_ll_to_xy(proj, grid2->sproj);

	int nhp = gcm->nhp(-1);
	// TODO: Consider copying iweightEp to a hash_map<int,double> for greater speed in loop.

	for (auto cell = gcm->grid1->cells_begin(); cell != gcm->grid1->cells_end(); ++cell) {
		if (gcm->masked1(cell->index)) continue;

		for (ihp=0; ihp<nhp; ++ihp) {

			double native_area = cell->area;
			double proj_area = area_of_proj_polygon(*cell, proj);

			long iE = cell->index;
			long iEp = iEp_base + iE;


			if (direction == ProjCorrect::NATIVE_TO_PROJ) {
				double coeff = native_area / proj_area;
				M.add(iEp, iE, coeff);
				// No weight vector needed here because it would be trivial.
				// In theory, the would match what's below, but things cancel out:
				// double coeff = weightE(iE) * native_area / proj_area;
				// M.add(iEp, iE, coeff);
				// weightE.add(iEp, weightE(iE));    // Only one thing gets added here, so no need for weight vector!

			} else {
				double coeff = iweightEp(iEp) * proj_area / native_area;
				M.add(iE, iEp, coeff);
				// After dividing by weightE, mean of final matrix will be (proj_area / native_area)
				weightE.add(iE, iweightEp(iEp));
			}
		}
	}
}





