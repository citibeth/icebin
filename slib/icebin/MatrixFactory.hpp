
namespace icebin {


struct RegridMatrix {
	CooMatrix M;		// B<-A; transpose for A<-B
	CooVector weightA, weightB;
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
Invert<TypeT> invert(TypeT const *M, bool inv)
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


template<TypeT>
class VectorAllocator {
	std::vector<std::unique_ptr<TypeT>> mem;
public:
	TypeT *alloc() {
		std::unique_ptr<TypeT> M(new TypeT());
		TypeT *ret = M.get();
		mem_matrix.push_back(std::move(M));
		return ret;
	}
}



class UrMatrices {
	VectorAllocator<SparseMatrix> matrix_mem;
	VectorAllocator<SparseVector> vector_mem;

	std::map<std::string, std::pair<SparseMatrix const *, bool>> regrids;
	std::map<std::string, std::pair<SparseVector const *, bool>> diags;


	/** Adds a regrid matrix and its variants.

	@param G The destination vector space of the matrix.  This must always be "G".
	@param X The source vectors space of the matrix.
	@param m_GvX Memory where to store the underlying matrix. */
	void add_regrids(std::string const &G, std::string const &X,
		std::function<void(SparseMatrix &)> const &regrid_fn)
	{
		// Make sure destination space is G
		if (G != "G") (*icebin_error)(-1,
			"Destination vector space for add_GvX must be \"G\" (exchnage grid)!");

		// Get the main matrix
		Sparse
		std::unique_ptr<SparseMatrix> M(new SparseMatrix());
		regrid_fn(*M);

		/* Since all raw matrices going into this have a destination G
		(exchange grid), they must all be sorted row-major.  The
		transpose matrices all go FROM G, so they must be sorted
		column-major.  But sorting the transpose of a matrix
		column-major is the same as sorting the original matrix row
		major.  So... all matrices must be sorted the same way
		(row-major) */
		M->consolidate({0,1});

		// Create the weight matrices
		std::unique_ptr<SparseVector> weightG(new SparseVector());
		std::unique_ptr<SparseVector> weightX(new SparseVector());
		for (auto ii=M->begin(); ii != M->end(); ++ii) {
			weightG->add({cell.index(0)}, cell.val());
			weightX->add({cell.index(1)}, cell.val());
		}
		weightG->consolidate({0});
		weightX->consolidate({0});

		std::string GvX = G+"v"+X;
		std::string XvG = X+"v"+G;

		regrids[GvX] = invert(&*M, false);
		regrids[XvG] = invert(&*M, true);			// true --> transpose
		diags["w"+GvX] = invert(&*weightG, false);
		diags["w"+XvG] = invert(&*weightX, false);
		diags["s"+GvX] = invert(&*weightG, true);	// true --> divide by this
		diags["s"+XvG] = invert(&*weightX, true);

		// Keep these matrices around before we go out of scope
		mem_matrix.push_back(std::move(M));
		mem_matrix.push_back(std::move(weightG));
		mem_matrix.push_back(std::move(weightX));
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

	void IceSheetUrMatrices(IceRegridder &sheet)
	{
		// Set up original matrices
		add_regrids("G", "Ep", std::bind(&IceRegridder::GvEp_noweight, _1));
		add_regrids("G", "Ap", std::bind(&IceRegridder::GvAp_noweight, _1));
		add_regrids("G", "I", std::bind(&IceRegridder::GvI_noweight, _1));

		SparseVector wEvEp(...);
		add_weight("E", "X", wEvX = product_ele(wEvEp, ur.diags["wEpvG"]);

		SparseVector wAvAp(...);
		add_weight("A", "X", wAvX = product_ele(wAvAp, ur.diags["wApvG"]);
	}
		
};

SparseMatrix prod(
	Invert<SparseVector> scalei,
	Invert<SparseMatrix> A,
	Invert<SparseVector> scalej,
	Invert<SparseMatrix> B,
	Invert<SparseVector> scalek)
{
	SparseMatrix ret;
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

Invert<SparseVector> NO_VECTOR

class IceSheetMatrices {
	IceSheetUrMatrices ur;

	SparseMatrix EvI, IvE;
	SparseMatrix AvI, IvA;
	SparseMatrix EvA, AvE;


	IceSheetMatrices(IceSheetUrMatrices &ur)
	{
		EvI = prod(ur.diags(
	}


	SparseMatrix AvI_noweight() {
		return prod(
			NO_VECTOR,
			ur.regrids["ApvG"],  ur.diags["sGvI"], ur.regrids["GvI"],
			NO_VECTOR);
	}

	SparseMatrix EvI_noweight() {
		return prod(
			NO_VECTOR,
			ur.regrids["EpvG"],  ur.diags["sGvI"], ur.regrids["GvI"],
			NO_VECTOR);
	}

	SparseMatrix IvA() {
		return prod(
			ur.diags["sIvG"],
			ur.regrids["IvG"], ur.diags["sGvAp"], ur.regrids["GvAp"],
			ur.diags["sApvA"]);
	}

	SparseMatrix IvE() {
		return prod(
			ur.diags["sIvG"],
			ur.regrids["IvG"], ur.diags["sGvEp"], ur.regrids["GvEp"],
			ur.diags["sEpvE"]);
	}


	/** Add to global inter-sheet weight vector, to be used for
	anything ending in A grid. */
	void wAvX(SparseVector &partial_cell_w)
		{ copy(partial_cell_w, ur.diags["wAvX"]); }
	void wEvX(SparseVector &partial_cell_w)
		{ copy(partial_cell_w, ur.diags["wEvX"]); }

	SparseMatrix EvA_noweight(SparseMatrix &EvA_global)
	{
		multiply(EvA_global,
			NO_VECTOR,
			ur.regrids["ApvG"], ur.diags["sGvEp"], ur.regrids["GvEp"],
			ur.diags["sEpvE"]);
	}

	SparseMatrix AvE_noweight(SparseMatrix &AvE_global)
	{
		multiply(AvE_global,
			NO_VECTOR,
			ur.regrids["EpvG"], ur.diags["sGvAp"], ur.regrids["GvAp"],
			ur.diags["sApvE"]);
	}

};

