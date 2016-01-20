enum class WeightType {
	NONE,			// Do not divide by area; result is in [kg] instead of [kg m-2]

	/** Divide by area of just portions of the cell overlapping with
	other relevant grids (eg: ice-covered portions).  Gives result in [kg m-2] */
	PARTIAL_CELL,

	WHOLE_CELL		// Divide by area of the ENTIRE cell, giving result in [kg m-2]
}

typedef spsparse::VectorCooArray<long, double, 2> CooMatrix;
typedef spsparse::VectorCooArray<long, double, 1> CooVector;

struct RegridMatrix {
	std::string A, B;	// String tags for A and B vector spaces.
	CooMatrix M;		// B<-A; transpose for A<-B
}

struct WeightVector {
	std::string A, B;	// String tags for A and B vector spaces.
	CooVector V;		// B<-A; invert for A<-B
}



template<class TypeT>
TypeT const &identity(TypeT *val) { return val; }


template<class TypeT>
class AltPtr {
	TypeT * urval;
	TypeT newval;
	TypeT *ptr;
public:
	AltPtr() : urval(0), ptr(0) {}
	AltPtr(TypeT *val) : urval(val), ptr(val) {}
	void set(TypeT &&_newval) {
		newval = _newval;
		ptr = &newval;
	}
	bool isset() { return ptr == &newval; }

	TypeT &operator*() { return *ptr; }
	TypeT *operator->() { return ptr; }
}


// -----------------------------------------------------
template<class TypeT>
class Invert {
	AltPtr<TypeT> M;
	bool invert;

	Invert() : M(0), invert(false) {}	// Null pointer
	Invert(TypeT const *_M, bool _invert);
	TypeT &operator*() { return M.operator*(); }
	TypeT *operator->() { return M.operator->(); }
};
template<class TypeT>
Invert<TypeT> invert(TypeT const &M, bool inv)
	{ return Invertible<TypeT>(M, inv); }

template<class TypeT>
Invert<TypeT>::Invert(TypeT const *_M, bool _invert) : M(_M), invert(_invert) {}

template<>
Invert<CooVector>::Invert(CooVector const *_M, bool _invert) : M(_M), invert(_invert)
{
	CooVector inv;
	for (auto ii=M->begin(); ii != M->end(); ++ii) inv.add(ii.indices(), 1.0/ii.val());
	M.set(std::move(inv));
}

// -----------------------------------------------------


template<class TypeT>
Invert<TypeT>::operator*() { return AltPtr<TypeT>(M); }

template<>
Invert<CooVector>::operator*() {
	AltPtr<TypeT> ptr;
	return AltPtr<TypeT>(M);
}




typedef Invert<CooVector> InvertVector;
class InvertMatrix : 


typedef std::pair<CooMatrix const *, bool> MatrixPtr;


class IceSheetUrMatrices {
	std::map<std::string, std::pair<CooMatrix const *, bool>> regrids;
	std::map<std::string, std::pair<CooVector const *, bool>> diags;

	RegridMatrix GvEp;		// Exchange <-- Elevation (Projected)
	RegridMatrix GvAp;		// Exchange <-- Atmosphere (Projected)
	RegridMatrix GvI;		// Exchange <-- Ice
	WeightMatrix wEvX;		// Inverse of scale matrix
	WeightMatrix wAvX;


	void add_regrids(RegridMatrix &BvA)
	{
		std::string BvA = B+"v"+A;
		std::string AvB = A+"v"+B;

		regrids[BvA] = invert(BvA.M, false);
		regrids[AvB] = invert(BvA.M, true);			// true --> transpose
		diags["w"+BvA] = invert(BvA.weightB, false);
		diags["w"+AvB] = invert(BvA.weightA, false);
		diags["s"+BvA] = invert(BvA.weightB, true);	// true --> divide by this
		diags["s"+AvB] = invert(BvA.weightA, true);
	}

	void add_weight(WeightVector &scale)
	{
		std::string BvA = B+"v"+A;
		std::string AvB = A+"v"+B;
		diags["s"+BvA] = invert(scale.V, false);
		diags["s"+AvB] = invert(scale.V, true);
		diags["w"+BvA] = invert(scale.V, true);
		diags["w"+AvB] = invert(scale.V, false);

	}

	void IceSheetUrMatrices(IceSheet &sheet)
	{
		// Set up original matrices
		add_regrids(GvEp = sheet->projelev_xx_iceexch(IceExch::EXCH));
		add_regrids(GvAp = sheet->projatm_xx_iceexch(IceExch::EXCH));
		add_regrids(GvI = sheet->GvI(IceExch::EXCH));

		CooVector wEvEp(...);
		add_weight("E", "X", wEvX = product_ele(wEvEp, ur.diags["wEpvG"]);

		CooVector wAvAp(...);
		add_weight("A", "X", wAvX = product_ele(wAvAp, ur.diags["wApvG"]);
	}
		
};

CooMatrix prod(
	Invert<CooVector> scalei,
	Invert<CooMatrix> A,
	Invert<CooVector> scalej,
	Invert<CooMatrix> B,
	Invert<CooVector> scalek)
{
	CooMatrix ret;
	spsparse::multiply(ret, 1.0,
		&*scalei,
		*A, A.invert,
		&*scalej,
		*B, B.invert,
		&*scalek,
		DuplicatePolicy::ADD,
		zero_nan);
	return ret;
}

Invert<CooVector> NO_VECTOR

class IceSheetMatrices {
	IceSheetUrMatrices ur;

	CooMatrix EvI, IvE;
	CooMatrix AvI, IvA;
	CooMatrix EvA, AvE;


	IceSheetMatrices(IceSheetUrMatrices &ur)
	{
		EvI = prod(ur.diags(
	}


	CooMatrix AvI() {
		return prod(
			NO_VECTOR,
			ur.regrids["ApvG"],  ur.diags["sGvI"], ur.regrids["GvI"],
			NO_VECTOR);
	}

	CooMatrix EvI() {
		return prod(
			NO_VECTOR,
			ur.regrids["EpvG"],  ur.diags["sGvI"], ur.regrids["GvI"],
			NO_VECTOR);
	}

	CooMatrix IvA_noweight() {
		return prod(
			ur.diags["sIvG"],
			ur.regrids["IvG"], ur.diags["sGvAp"], ur.regrids["GvAp"],
			ur.diags["sApvA"]);
	}

	CooMatrix IvE_noweight() {
		return prod(
			ur.diags["sIvG"],
			ur.regrids["IvG"], ur.diags["sGvEp"], ur.regrids["GvEp"],
			ur.diags["sEpvE"]);
	}


	/** Add to global inter-sheet weight vector, to be used for
	anything ending in A grid. */
	void wAvX(CooVector &partial_cell_w)
		{ copy(partial_cell_w, ur.diags["wAvX"]); }
	void wEvX(CooVector &partial_cell_w)
		{ copy(partial_cell_w, ur.diags["wEvX"]); }

	CooMatrix EvA_noweight(CooMatrix &EvA_global)
	{
		multiply(EvA_global,
			ur.diags["sAwX"],
			ur.regrids["ApvG"], ur.diags["sGvEp"], ur.regrids["GvEp"],
			ur.diags["sEpvE"]);
	}

	CooMatrix AvE_noweight(CooMatrix &AvE_global)
	{
		multiply(AvE_global,
			ur.diags["sEwX"],
			ur.regrids["EpvG"], ur.diags["sGvAp"], ur.regrids["GvAp"],
			ur.diags["sApvE"]);
	}

};

