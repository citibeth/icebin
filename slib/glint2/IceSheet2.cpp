



/** Adds to a sparse matrix that projects from atmosphere to projected
atmosphere.  Note there is one projected atmosphere per ice sheet.

If NATIVE_TO_PROJ (A --> Ap):
    fAp = M fA    for each ice sheet
If PROJ_TO_NATIVE (Ap --> A):
	fA = diag(1/weightA) sum_{ice sheets}(M fAp)
*/
void atm_project(			// Transformation: Ap<-A
giss::CooVector &A2Ap,
giss::CooVecotr &Ap2A,
giss::MapSparseVector<int,double> const *iweightAp,		// IN : Weight vector from iceexch_xx_projatm: portion of this projected atm cell used
giss::MapSparseVector<int,double> weightA)				// OUT: Scaling to use when applying this transformation
{
	giss::Proj2 proj;
	gcm->grid1->get_ll_to_xy(proj, grid2->sproj);

	// TODO: Consider copying iweightAp to a hash_map<int,double> for greater speed in loop.
	for (auto cell = gcm->grid1->cells_begin(); cell != gcm->grid1->cells_end(); ++cell) {
		if (gcm->masked1(cell->index)) continue;

		double native_area = cell->area;
		double proj_area = area_of_proj_polygon(*cell, proj);

		index_t iA = cell->index;

		// Ap<-A
		double coeff = native_area / proj_area;
		A2Ap.add(iAp, iA, coeff);

		// A<-Ap
		coeff = iweightAp(iAp) * proj_area / native_area;
		Ap2A.add(iA, iAp, coeff);
		// After dividing by weightA, mean of final matrix will be (proj_area / native_area)
		weightA.add(iA, iweightAp(iAp));
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

			long iE = hc_index->ik_to_index(cell->index, ihp);
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


/** Parse a matrix spec */
class MatrixSpec {

	std::string src, dst;
	array<int, 2> sort_order;	

	MatrixSpec(std::string const &sspec) : sort_order({{-1, -1}})
	{
		int arrow = sspec.find("<-");
		if (arrow < 0) {
			src = sspec;
			return;
		}

		size_t slen = sspec.size()-1;
		src = sspec.substr(0,arrow);
		if (sspec[slen] == '|') {
			dst = sspec.substr(arrow+2, (slen-1)-(arrow+2));
			sort_order = COLUMN_MAJOR;
		} else if (sspec[sspec.size()-1] == '-') {
			dst = sspec.substr(arrow+2, (slen-1)-(arrow+2));
			sort_order = ROW_MAJOR;
		} else {
			dst = sspec.substr(arrow+2, slen-(arrow+2));
		}
	}

	/** Returns 1 if the two specs are the same (eg: A<-B), and -1 if they
	are reversed (B<-A).  Throws an error if they don't match at all. */
	int direction(MatrixSpec &other)
	{
		if ((src == other.src) && (dst == other.dst)) return 1;
		if ((src == other.dst) && (dst == other.src)) return -1;
		fprintf(stderr, "Mismatched matrix spec: (%s<-%s) vs (%s<-%s)\n",
			other.src, other.dst, src, dst);
		giss::error(-1);
	}
}






/** A matrix combined with a scale vector: M: A -> B, M=(M0, w)
    The definition is:
        Mx = diag(scale) M0x
NOTES:

 1. The scale vector can be considered to be the inverse of the area
    contribution of all grid cells in A to each grid cell in B.

 2. The scale vector might be shared among many MatrixAndScale.
    For example... M might regrid from an individual ice sheet to a
    global GCM grid.  The scale would be based on the contributions
    of ALL the ice sheets to the global grid.
*/
template<class IndexT, class ValueT>
class ScaledMatrix {
	CooMatrix<IndexT, ValueT> M;
	CooVector<IndexT, ValueT> scale;
};

class RegridMatrix {
	MatrixSpec spec;
	glint2::CooMatrix A2B;		// A2B   (B,A)
	std::unique_ptr<glint2::CooMatrix> &B2A;
	glint2::CooVector weightA;		// B
	glint2::CooVector weightB;		// A

public:
	RegridMatrix(std::string sspec) : spec(MatrixSpec(sspec))
		{}

	/** Remove derived stuff */
	void clean()
		{ B2A.reset(); }

	glint2::CooMatrix &operator()(int direction = 1,
		std::array<std::int, 2> sort_order = {{-1,-1}})
	{
		if (direction > 0) {
			// Return the forward matrix
			if (sort_order[0] >= 0 && sort_order != A2B.sort_order)
				A2B.consolidate(sort_order);
			return A2B;
		} else {
			if (!B2A.get()) B2A.reset(new glint2::CooMatrix(transpose(A2B)));
			if (sort_order[0] >= 0 && sort_order != B2A->sort_order)
				B2A->consolidate(sort_order);
			return *B2A;
		}
	}

	void glint2::CooMatrix &operator()(std::string const &sspec)
	{
		ospec = MatrixSpec(sspec);
		return (*this)(ospec.direction(spec), ospec.sort_order);
	}

	glint2::CooVector &weight(std::string const &sgrid)
	{
		if (sgrid == spec.src) return weightA;
		if (sgrid == spec.dst) return weightB;
		fprintf(stderr, "Grid spec %s doesn't match anything in (%s<-%s)\n", sgrid, spec.dst, spec.src);
		giss::exit(-1);
	}

	glint2::CooVector scale(std::string const &sgrid)
	{
		glint2::CooVector &w(weight(sgrid));
		return inverse_ele(w);
	}
};


class ProjMatrix {
	CooVector<index_t, double> A2Ap;		// diag(A2Ap) is the matrix
	CooVector<index_t, double> Ap2A;		// diag(Ap2A) is the matrix == diag(1/A2Ap)

	CooVector<index_t, double> weightA;
	// CooVector<index_t, double> weightAp;		// Always 1.0

public:
	ProjMatrix(std::string sspec) : spec(MatrixSpec(sspec))
		{}

	void glint2::CooMatrix &operator()(int direction = 1)
	{
		if (direction < 0) return Ap2A;
		else return A2Ap;
	}

	void glint2::CooMatrix &operator()(std::string const &sspec)
	{
		ospec = MatrixSpec(sspec);
		return (*this)(ospec.direction(spec), ospec.sort_order);
	}


	glint2::CooVector &weight(std::string const &sgrid)
	{
		if (sgrid == spec.src) return weightA;
		fprintf(stderr, "Grid spec %s doesn't match anything in (%s<-%s)\n", sgrid, spec.dst, spec.src);
		giss::exit(-1);
	}

	glint2::CooVector scale(std::string const &sgrid)
	{
		glint2::CooVector &w(weight(sgrid));
		return inverse_ele(w);
	}


};


void IceSheet::regrid_matrix(std::string const &spec)
{
}


template<
	template<class KeyTT, class ValTT> class BaseTpl,
	class KeyT, class ValT>
class UniquePtrDict : public BaseTpl<KeyT, std::unique_ptr<ValT>>
{
	void inset(typename MapT::key_type const &key, typename MapT::mapped_type &&val)
		{ super::insert(key, std::unique_ptr<ValT>(new ValT(std::move(val)))); }

	void inset(typename MapT::key_type const &key, typename MapT::mapped_type &val)
		{ super::insert(key, std::unique_ptr<ValT>(new ValT(val))); }

	void inset(typename MapT::key_type const &key, typename MapT::mapped_type val)
		{ super::insert(key, std::unique_ptr<ValT>(new ValT(std::move(val)))); }

	void inset(typename MapT::key_type const &key, std::unique_ptr<ValT> &&val)
		{ super::insert(key, std::move(val); }
};


class CouplingMatrices
{
public:

	UniquePtrDict<std::map, std::string, giss::CooMatrix> regrids;
	UniquePtrDict<std::map, std::string, giss::CooMatrix> weights;

	giss::CooMatrix &operator()(std::string const &sspec)
		{ return regrids.at(sspec); }

	giss::CooMatrix &weight()(std::string const &sspec)
		{ return weights.at(sspec); }

	void CouplingMatrices(IceSheet *sheet)
	{
		// --------------- Base Matrices
		// Base regrid matrices
		RegridMatrix M("G<-Ep");	// Ep->G
		sheet->projelev_xx_iceexch(M("G<-Ep"), M.weight("Ep"), M.weight("G"), IceExch::EXCH);

		RegridMatrix R("G<-Ap");	// Ap->G
		sheet->iceexch_xx_projatm(R("Ap<-G"), R.weight("Ap"), R.weight("G"), IceExch::EXCH);
		RegridMatrix X("G<-I");
		sheet->ice_xx_exch(X("G<-I"), X.weight("I"), R.weight("G"));

		// Base projection matrices
		ProjMatrix ProjE("Ep<-E");
		sheet->atm_project(ProjE("Ep<-E"), ProjE("E<-Ep"), M.weight("Ep"), ProjE.weight("E"));
		ProjMatrix ProjA("Ap<-A");
		sheet->atm_project(ProjA("Ap<-A"), ProjA("A<-Ap"), R.weight("Ap"), ProjA.weight("A"));


		// ---------- Derived matrices
		// ---- E2I
		CooVector<index_t, double> E_scalei = prod_ele(ProjE("E<-Ep"), M.scale("Ep"));
		regrids.insert("I<-E",
			prod(&X.scale("I"), X("I<-G-"), &M.scale("G"), M("G<-Ep|"), &ProjE("Ep<-E")));

		// ---- I2E (requires post-application of E_weight)
		regrids.insert("E<-I",
			prod(&E_scalei, M("Ep<-G-"), &X.scale("G"), X("G<-I|"), NULL);


		// ---- A2I  (symmetric with E2I)
		CooVector<index_t, double> A_scalei = prod_ele(ProjA("A<-Ap"), R.scale("Ap"));
		regrids.insert("I<-A",
			prod(A2I, &X.scale("I"), X("I<-G-"), &R.scale("G"), R("G<-Ap|"), &ProjA("Ap<-A"));

		// ---- I2A (requires post-application of A_weight)
		regrids.insert("A<-I"
			prod(A_scalei, R("Ap<-G-"), &X.scale("G"), X("G<-I|"), NULL);



		// ---- A2E (requires post-application of E_weight)
		regrids.insert("E<-A",
			prod(&E_scalei, M("Ep<-G-"), &X.scale("G"), R("G<-Ap|"), &ProjA("Ap<-A"));

		// ---- E2A (requires post-application of E_weight)
		regrids.insert("A<-E",
			prod(&A_scalei, R("Ap<-G-"), &X.scale("G"), M("G<-Ep|"), &ProjE("Ep<-E"));


		// Remove our weight matrices from temporaries
		weights.insert("A", std::move(ProjA.weight("A")));
		weights.insert("E", std::move(ProjA.weight("A")));
	}
		
