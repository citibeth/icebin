/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <giss/CooVector.hpp>
#include <giss/ncutil.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <giss/IndexTranslator.hpp>
#include <giss/IndexTranslator2.hpp>
#include <galahad/qpt_c.hpp>
#include <galahad/eqp_c.hpp>
#include <giss/ncutil.hpp>
#include <giss/enum.hpp>

namespace glint2 {

void MatrixMaker::clear()
{
	sheets.clear();
	sheets_by_id.clear();
	grid1.reset();
	mask1.reset();
	hpdefs.clear();
	// hcmax.clear();
}

void MatrixMaker::realize() {

	// ---------- Check array bounds
	long n1 = grid1->ndata();
	if (mask1.get() && mask1->extent(0) != n1) {
		fprintf(stderr, "mask1 for %s has wrong size: %d (vs %d expected)\n",
			mask1->extent(0), n1);
		throw std::exception();
	}

	// ------------- Realize the ice sheets
	for (auto sheet=sheets.begin(); sheet != sheets.end(); ++sheet)
		sheet->realize();

	// ------------- Set up HCIndex
	hc_index = HCIndex::new_HCIndex(_hc_index_type, *this);
}

int MatrixMaker::add_ice_sheet(std::unique_ptr<IceSheet> &&sheet)
{
	if (sheet->name == "") {
		fprintf(stderr, "MatrixMaker::add_ice_sheet(): Sheet must have a name\n");
		throw std::exception();
	}

	int const index = _next_sheet_index++;
	sheet->index = index;
printf("MatrixMaker: %p.sheetno = %d\n", &*sheet, sheet->index);
	sheet->gcm = this;
	
	sheets_by_id.insert(std::make_pair(sheet->index, sheet.get()));
	sheets.insert(sheet->name, std::move(sheet));
	return index;
}

// --------------------------------------------------------------
/** NOTE: Allows for multiple ice sheets overlapping the same grid cell (as long as they do not overlap each other, which would make no physical sense). */
void MatrixMaker::fgice(giss::CooVector<int,double> &fgice1)
{

	// Accumulate areas over all ice sheets
	giss::SparseAccumulator<int,double> area1_m_hc;
	fgice1.clear();
	for (auto sheet = sheets.begin(); sheet != sheets.end(); ++sheet) {

		// Local area1_m just for this ice sheet
		giss::SparseAccumulator<int,double> area1_m;
		sheet->accum_areas(area1_m);

		// Use the local area1_m to contribute to fgice1
		giss::Proj2 proj;
		grid1->get_ll_to_xy(proj, sheet->grid2->sproj);
		for (auto ii = area1_m.begin(); ii != area1_m.end(); ++ii) {
			int const i1 = ii->first;
			double ice_covered_area = ii->second;
			Cell *cell = grid1->get_cell(i1);
			if (!cell) continue;	// Ignore cells in the halo
			double area1 = area_of_proj_polygon(*cell, proj);
			fgice1.add(i1, ice_covered_area / area1);

		}
	}
	fgice1.sort();
}

// --------------------------------------------------------------
/** Change this to a boost function later
@param trans_2_2p Tells us which columns of conserv (i2) are active,
       after masking with mask1, mask1h and mask2.
*/
static std::unique_ptr<giss::VectorSparseMatrix> remove_small_constraints(
giss::VectorSparseMatrix const &in_constraints_const,
int min_row_count)
{
	// Const cast because we don't know how to do const_iterator() right in VectorSparseMatrix
	auto in_constraints(const_cast<giss::VectorSparseMatrix &>(in_constraints_const));

	std::set<int> delete_row;		// Rows to delete
	std::set<int> delete_col;		// Cols to delete

	// Make sure there are no constraints (rows) with too few variables (columns).
	// Uses an iterative process
	std::vector<int> row_count(in_constraints.nrow);
	for (;;) {
		// Count rows
		row_count.clear(); row_count.resize(in_constraints.nrow);
		for (auto oi = in_constraints.begin(); oi != in_constraints.end(); ++oi) {
			int i2 = oi.col();

			// Loop if it's already in our delete_row and delete_col sets
			if (delete_row.find(oi.row()) != delete_row.end()) continue;
			if (delete_col.find(i2) != delete_col.end()) continue;

			++row_count[oi.row()];
		}


		// Add to our deletion set
		int num_deleted = 0;
		for (auto oi = in_constraints.begin(); oi != in_constraints.end(); ++oi) {
			int i2 = oi.col();

			// Loop if it's already in our delete_row and delete_col sets
			if (delete_row.find(oi.row()) != delete_row.end()) continue;
			if (delete_col.find(i2) != delete_col.end()) continue;

			if (row_count[oi.row()] < min_row_count) {
				++num_deleted;
				delete_row.insert(oi.row());
				delete_col.insert(i2);
			}
		}

		// Terminate if we didn't remove anything on this round
printf("num_deleted = %d\n", num_deleted);
		if (num_deleted == 0) break;
	}


	// Copy over the matrix, deleting rows and columns as planned
	std::unique_ptr<giss::VectorSparseMatrix> out_constraints(
		new giss::VectorSparseMatrix(giss::SparseDescr(in_constraints)));
	for (auto oi = in_constraints.begin(); oi != in_constraints.end(); ++oi) {
		int i2 = oi.col();

		// Loop if it's already in our delete_row and delete_col sets
		if (delete_row.find(oi.row()) != delete_row.end()) continue;
		if (delete_col.find(i2) != delete_col.end()) continue;

		out_constraints->set(oi.row(), i2, oi.val());
	}
	return out_constraints;
}
// -------------------------------------------------------------
/** Checksums an interpolation matrix, to ensure that the sume of weights
for each output grid cell is 1.  This should always be the case, no
matter what kind of interpolation is used. */
static bool checksum_interp(giss::VectorSparseMatrix &mat, std::string const &name,
double epsilon)
{
	bool ret = true;
	auto rowsums(mat.sum_per_row_map());
	for (auto ii = rowsums.begin(); ii != rowsums.end(); ++ii) {
		int row = ii->first;
		double sum = ii->second;
		if (std::abs(sum - 1.0d) > epsilon) {
			printf("rowsum != 1 at %s: %d %g\n", name.c_str(), row, sum);
			ret = false;
		}
	}
	return false;
}
// -------------------------------------------------------------
class I3XTranslator {
	HCIndex *hc_index;
	int nhc;

public:

	void init(HCIndex *_hc_index, int _nhc)
	{
		hc_index = _hc_index;
		nhc = _nhc;
	}

	// Identity Transformation
	void init_identity() {
		hc_index = 0;
		nhc = -1;
	}

	I3XTranslator() : hc_index(0), nhc(-1) {}

	int i3x_to_i3(int i3x)
	{
		if (!hc_index) return i3x;		// Identity

		int i1 = i3x / nhc;
		int k = i3x - i1 * nhc;
		return hc_index->ik_to_index(i1, k);
	}

	int i3_to_i3x(int i3)
	{
		if (!hc_index) return i3;		// Identity

		int i1, k;
		hc_index->index_to_ik(i3, i1, k);
		int i3x = i1 * nhc + k;
	//printf("i3=%d (%d, %d) --> i3x=%d\n", i3, i1, k, i3x);
		return i3x;
	}
};
// -------------------------------------------------------------
struct UsedAll {
	std::set<int> used1;
	std::set<std::pair<int,int>> used4;
	std::set<int> used3;
	std::set<int> used3x;
	giss::IndexTranslator trans_1_1p;
	giss::IndexTranslator2 trans_4_4p;
	giss::IndexTranslator trans_3x_3p;

	std::unique_ptr<giss::VectorSparseMatrix> RMp;
	std::unique_ptr<giss::VectorSparseMatrix> Sp;
	std::unique_ptr<giss::VectorSparseMatrix> XMp;
	blitz::Array<double,1> area1p_inv;

	UsedAll() :
		trans_1_1p("trans_1_1p"),
		trans_4_4p("trans_4_4p"),
		trans_3x_3p("trans_3x_3p") {}

};

class GetSubID {
	QPAlgorithm _qp_algorithm;
public:
	GetSubID(QPAlgorithm qp_algorithm) : _qp_algorithm(qp_algorithm) {}
	int operator()(int i1) {
		if (_qp_algorithm == QPAlgorithm::SINGLE_QP) return 0;
		return i1;
	}
};

static int get_subid(int i1) { return i1; }



/** @params f2 Some field on each ice grid (referenced by ID).  Do not have to be complete. */
giss::CooVector<int, double>
MatrixMaker::iceinterp_to_hp(
std::map<int, blitz::Array<double,1>> &f2_or_4s,		// Actually f2 or f4
blitz::Array<double,1> initial3,
IceInterp src,
QPAlgorithm qp_algorithm)
{
printf("BEGIN MatrixMaker::iceinterp_to_hp()\n");

	// Convert our input from ice grid to interpolation grid
	// (but ONLY if needed)
	std::map<int, blitz::Array<double,1>> _f4s_new;
	std::map<int, blitz::Array<double,1>> *f4s;
	if (src == IceInterp::INTERP) {
		// f2_or_4s is already in interpolation grid space
		f4s = &f2_or_4s;
	} else {
		f4s = &_f4s_new;
		for (auto f2i=f2_or_4s.begin(); f2i != f2_or_4s.end(); ++f2i) {
			IceSheet *sheet = (*this)[f2i->first];
			blitz::Array<double,1> &f2(f2i->second);
			f4s->insert(std::make_pair(f2i->first, sheet->ice_to_interp(f2)));
		}
	}


	bool const convert_3_3x = false;		// Didn't seem to help

	// =============== Partition big QP problems into many little ones
	// (if caller requested)
	GetSubID get_subid(qp_algorithm);
	// Temporary variables required on a per-sub-problem basis
	giss::MapDict_Create<int, UsedAll> used;		// subid -> A bunch of used sets

	// =============== Set up basic vector spaces for optimization problem

	// Used in constraints
	std::unique_ptr<giss::VectorSparseMatrix> RM0(hp_to_atm());	// 3->1
	for (auto ii = RM0->begin(); ii != RM0->end(); ++ii) {
		int i1 = ii.row();
		UsedAll &ua(*used[get_subid(i1)]);

		// Check RM is local for MULTI_QP algorithm
		int i3 = ii.col();
		if (qp_algorithm == QPAlgorithm::MULTI_QP) {
			int i1b, k;
			hc_index->index_to_ik(i3, i1b, k);
			if (i1b != i1) {
				fprintf(stderr, "RM (hp2atm) matrix is non-local!\n");
				throw std::exception();
			}
		}

		ua.used1.insert(i1);
		ua.used3.insert(i3);
	}

// In some cases in the past, QP optimization has not worked well
// when there are grid cells with very few entries in the
// constraints matrix.  Not an issue here now.
#if 0
	std::unique_ptr<giss::VectorSparseMatrix> RM(
		remove_small_constraints(*RM0, 2));
	RM0.reset();
#else
	std::unique_ptr<giss::VectorSparseMatrix> RM(std::move(RM0));
#endif

	giss::SparseAccumulator<int,double> area1;
	giss::MapDict<int, giss::VectorSparseMatrix> Ss;
	giss::MapDict<int, giss::VectorSparseMatrix> XMs;
	std::map<int, size_t> size4;	// Size of each ice vector space
	for (auto f4i=f4s->begin(); f4i != f4s->end(); ++f4i) {
		IceSheet *sheet = (*this)[f4i->first];

		std::unique_ptr<giss::VectorSparseMatrix> S(
			sheet->iceinterp_to_projatm(area1, IceInterp::INTERP));		// 4 -> 1
		if (correct_area1) S = multiply(
			*sheet->atm_proj_correct(ProjCorrect::PROJ_TO_NATIVE),
			*S);
		for (auto ii = S->begin(); ii != S->end(); ++ii) {
			int i1 = ii.row();
			UsedAll &ua(*used[get_subid(i1)]);
			ua.used1.insert(i1);
			ua.used4.insert(std::make_pair(sheet->index, ii.col()));
		}

		std::unique_ptr<giss::VectorSparseMatrix> XM(
			sheet->hp_to_iceinterp(IceInterp::INTERP));				// 3 -> 4
//		checksum_interp(*XM, "XM");

		for (auto ii = XM->begin(); ii != XM->end(); ++ii) {
			int i4 = ii.row();

			int i3 = ii.col();
			int i1, k;
			hc_index->index_to_ik(i3, i1, k);

			UsedAll &ua(*used[get_subid(i1)]);
			ua.used4.insert(std::make_pair(sheet->index, i4));
			ua.used3.insert(i3);
		}
printf("MatrixMaker::ice_to_hp() 4\n");

		size4[sheet->index] = sheet->n4();

		// Store away for later reference
		Ss.insert(sheet->index, std::move(S));
		XMs.insert(sheet->index, std::move(XM));
	}


	// -------- Set up the i3 <-> i3x transformation (if we want to)
	I3XTranslator trans_3_3x;
	if (convert_3_3x) {
		int max_k = 0;		// Maximum height class index
		for (auto ua = used.begin(); ua != used.end(); ++ua) {

			// Convert from i3 to i3x (renumbered height class indices)
			for (auto p3 = ua->used3.begin(); p3 != ua->used3.end(); ++p3) {
				int i1, k;
				hc_index->index_to_ik(*p3, i1, k);
				max_k = std::max(k, max_k);
			}
		}

		trans_3_3x.init(&*hc_index, max_k + 1);

		for (auto ua = used.begin(); ua != used.end(); ++ua) {

			for (auto p3 = ua->used3.begin(); p3 != ua->used3.end(); ++p3) {
				int i3x = trans_3_3x.i3_to_i3x(*p3);
				ua->used3x.insert(i3x);
			}
			ua->used3.clear();		// No longer needed
		}
	} else {
		trans_3_3x.init_identity();
		for (auto ua = used.begin(); ua != used.end(); ++ua) {
			ua->used3x = std::move(ua->used3);
		}
	}


	// -------------- Set up destination renumbered matrices
	for (auto ua = used.begin(); ua != used.end(); ++ua) {
		ua->trans_1_1p.init(n1(), ua->used1);
		ua->trans_4_4p.init(&size4, ua->used4);
		ua->trans_3x_3p.init(n3(), ua->used3x);

		int n1p = ua->trans_1_1p.nb();
		int n4p = ua->trans_4_4p.nb();
		int n3p = ua->trans_3x_3p.nb();

		// Translate to new matrices
		ua->RMp.reset(new giss::VectorSparseMatrix(giss::SparseDescr(n1p, n3p)));
		ua->Sp .reset(new giss::VectorSparseMatrix(giss::SparseDescr(n1p, n4p)));
		ua->XMp.reset(new giss::VectorSparseMatrix(giss::SparseDescr(n4p, n3p)));
		ua->area1p_inv.resize(n1p);
	}

	// ----------------- Copy data into those matrices
printf("Translating RM\n");
	for (auto ii = RM->begin(); ii != RM->end(); ++ii) {
		int i1 = ii.row();
		UsedAll *ua(used[get_subid(i1)]);

		int i3x = trans_3_3x.i3_to_i3x(ii.col());
		ua->RMp->add(
			ua->trans_1_1p.a2b(i1),
			ua->trans_3x_3p.a2b(i3x), ii.val());
	}


	for (auto f4i=f4s->begin(); f4i != f4s->end(); ++f4i) {
		int const index = f4i->first;
		IceSheet *sheet = (*this)[index];

		// Source matrices
		giss::VectorSparseMatrix *S(Ss[index]);
		giss::VectorSparseMatrix *XM(XMs[index]);

printf("Translating S: %d\n", index);
		for (auto ii = S->begin(); ii != S->end(); ++ii) {
			int i1 = ii.row();
			UsedAll *ua(used[get_subid(i1)]);
			ua->Sp->add(
				ua->trans_1_1p.a2b(i1),
				ua->trans_4_4p.a2b(std::make_pair(index, ii.col())),
				ii.val());
		}

printf("Translating XM: %d\n", index);
		for (auto ii = XM->begin(); ii != XM->end(); ++ii) {
			int i4 = ii.row();

			int i3 = ii.col();
			int i3x = trans_3_3x.i3_to_i3x(i3);
			int i1, k;
			hc_index->index_to_ik(i3, i1, k);
			UsedAll *ua(used[get_subid(i1)]);

			ua->XMp->add(
				ua->trans_4_4p.a2b(std::make_pair(index, ii.row())),
				ua->trans_3x_3p.a2b(i3x),
				ii.val());
		}
	}

	// ----------- Translate area1 -> area1p
	for (auto ua = used.begin(); ua != used.end(); ++ua) {
		int n1p = ua->trans_1_1p.nb();

		ua->area1p_inv.resize(n1p);
		ua->area1p_inv = 0;
	}

	for (auto ii = area1.begin(); ii != area1.end(); ++ii) {
		int i1 = ii->first;
		UsedAll *ua(used[get_subid(i1)]);
		int i1p = ua->trans_1_1p.a2b(i1);
		ua->area1p_inv(i1p) += ii->second;
	}

	for (auto ua = used.begin(); ua != used.end(); ++ua) {
		int n1p = ua->trans_1_1p.nb();

		for (int i1p=0; i1p<n1p; ++i1p) {
			if (ua->area1p_inv(i1p) != 0) ua->area1p_inv(i1p) = 1.0d / ua->area1p_inv(i1p);
		}

		// ---------- Divide Sp by area1p to complete the regridding matrix
		for (auto ii = ua->Sp->begin(); ii != ua->Sp->end(); ++ii) {
			int i1p = ii.row();
			ii.val() *= ua->area1p_inv(i1p);
		}
	}

	giss::CooVector<int, double> ret3;		// Function return value
	for (auto ua = used.begin(); ua != used.end(); ++ua) {
		int n1p = ua->trans_1_1p.nb();
		int n4p = ua->trans_4_4p.nb();
		int n3p = ua->trans_3x_3p.nb();
printf("--------------------- QP %d: n1p=%d, n4p=%d, n3p=%d\n", ua.key(), n1p, n4p, n3p);

		// -------- Translate f4 -> f4p
		// Ignore elements NOT listed in the translation
		blitz::Array<double,1> f4p(n4p);
		f4p = 0;
		for (int i4p = 0; i4p < n4p; ++i4p) {
			std::pair<int,int> const &a(ua->trans_4_4p.b2a(i4p));
			int index = a.first;
			int i4 = a.second;
			f4p(i4p) = (*f4s)[index](i4);
		}


	 	// ---------- Allocate the QPT problem
		// m = # constraints = n1p (size of atmosphere grid)
		// n = # variabeles = n3p
		galahad::qpt_problem_c qpt(n1p, n3p, true);

		// ================ Objective Function
		// 1/2 (XM F_E - F_I)^2    where XM = (Ice->Exch)(Elev->Ice)
		// qpt%H = (XM)^T (XM),    qpt%G = f_I \cdot (XM),        qpt%f = f_I \cdot f_I

		// -------- H = 2 * XMp^T XMp
		giss::VectorSparseMatrix XMp_T(giss::SparseDescr(ua->XMp->ncol, ua->XMp->nrow));
		transpose(*ua->XMp, XMp_T);
		std::unique_ptr<giss::VectorSparseMatrix> H(multiply(XMp_T, *ua->XMp));	// n3xn3

		// Count items in H lower triangle
		size_t ltri = 0;
		for (auto ii = H->begin(); ii != H->end(); ++ii)
			if (ii.row() >= ii.col()) ++ltri;

		// Copy ONLY the lower triangle items to GALAHAD
		// (otherwise, GALAHAD won't work)
printf("qpt.alloc_H(%d)\n", ltri);
		qpt.alloc_H(ltri);
		giss::ZD11SparseMatrix H_zd11(qpt.H, 0);
		for (auto ii = H->begin(); ii != H->end(); ++ii) {
			if (ii.row() >= ii.col()) {
				H_zd11.add(ii.row(), ii.col(), 2.0d * ii.val());
			}
		}

printf("AA1 Done\n");
		// -------- Linear term of obj function
		// G = -2*f4p \cdot XMp
		for (int i=0; i < qpt.n; ++i) qpt.G[i] = 0;
		for (auto ii = ua->XMp->begin(); ii != ua->XMp->end(); ++ii) {
			qpt.G[ii.col()] -= 2.0d * f4p(ii.row()) * ii.val();
		}

		// --------- Constant term of objective function
		// f = f4p \cdot f4p
		qpt.f = 0;
		for (int i4p=0; i4p<n4p; ++i4p) {
			qpt.f += f4p(i4p) * f4p(i4p);
		}

		// De-allocate...
//		H.reset();
//		ua->XMp->clear();
//		XMp_T.clear();

		// ============================ Constraints
		// RM x = Sp f4p

		// qpt.A = constraints matrix = RMp
		qpt.alloc_A(ua->RMp->size());
		giss::ZD11SparseMatrix A_zd11(qpt.A, 0);
		copy(*ua->RMp, A_zd11);

		// Constraints: Ax + C = 0
		// qpt.C = equality constraints RHS = Sp * f4p
		for (int i=0; i<n1p; ++i) qpt.C[i] = 0;
		for (auto ii = ua->Sp->begin(); ii != ua->Sp->end(); ++ii) {
			int i1p = ii.row();		// Atm
			int i4p = ii.col();		// Ice
			qpt.C[i1p] -= f4p(i4p) * ii.val();
		}

		// De-allocate
//		ua->RMp.reset();	// ->clear();

		// =========================== Initial guess at solution
		bool nanerr = false;
		bool zero_initial(initial3.size() == 0);
		for (int i3p=0; i3p<n3p; ++i3p) {
			if (zero_initial) {
				// No initial guess, start at 0
				qpt.X[i3p] = 0;
			} else {
				// Use initial guess supplied by user
				int i3x = ua->trans_3x_3p.b2a(i3p);
				int i3 = trans_3_3x.i3x_to_i3(i3x);
				double val = initial3(i3);
				if (std::isnan(val)) {
					fprintf(stderr, "ERROR: ice_to_hp(), NaN in initial guess, i3=%d\n", i3);
					nanerr = true;
					qpt.X[i3p] = 0;
				} else {
					qpt.X[i3p] = val;
				}
			}
		}

		// =========================== Solve the Problem!
		double infinity = 1e20;
		eqp_solve_simple(qpt.this_f, infinity);

		// ========================================================
		// ========================================================

		// --------- Pick out the answer and convert back to standard vector space
		for (int i3p=0; i3p<n3p; ++i3p) {
			int i3x = ua->trans_3x_3p.b2a(i3p);
			int i3 = trans_3_3x.i3x_to_i3(i3x);
			ret3.add(i3, qpt.X[i3p]);
		}
	}

	return ret3;
}


// --------------------------------------------------------------
/** @param force_lambda If true, then just return $f3 = \Lambda f1$,
           the physicaly most correct answer.
*/
giss::CooVector<int, double> MatrixMaker::atm_to_hp(blitz::Array<double,1> f1,
bool force_lambda)
{
	// RM = hp --> atm conversion
	std::unique_ptr<giss::VectorSparseMatrix> RM0(hp_to_atm());	// 3->1

	// Find non-zero rows and columns of RM.
	// Also see if RM is local (which eliminates need for quad opt)
	std::unordered_map<int,double> sum1;	// Sum of RM elements
	std::set<int> used1;
	std::set<int> used3;
	bool rm_local = true;
	for (auto ii = RM0->begin(); ii != RM0->end(); ++ii) {
		int i1 = ii.row();

		// Check RM is local
		int i3 = ii.col();

		if (rm_local) {
			int i1b, k;
			hc_index->index_to_ik(i3, i1b, k);
			if (i1b != i1) rm_local = false;
		}

		sum1[i1] += ii.val();
		used1.insert(i1);
		used3.insert(i3);
	}

	// Compute 1 / sum1, this is the factor to multiply by.
	// sum1(RM) = x means that height point = 1 gets scaled to x in atm grid.
	// So we must divide by x to rescale from atm to height point
	for (auto ii=sum1.begin(); ii!=sum1.end(); ++ii) {
printf("sum1[%d] = %g\n", ii->first, ii->second);
		ii->second = 1.0d / ii->second;
	}
	std::unordered_map<int,double> sum1_inv(std::move(sum1));

	// For local RM: Answer is \Lambda f1
	// (i.e. just repeat each atmosphere value for each elevation point)
printf("atm_to_hp: rm_local = %d\n", rm_local);
	if (force_lambda || rm_local) {
		giss::CooVector<int, double> ret;
		for (auto p3 = used3.begin(); p3 != used3.end(); ++p3) {
			int i3 = *p3;
			int i1, k;
			hc_index->index_to_ik(i3, i1, k);
			ret.add(i3, f1(i1) * sum1_inv[i1] );
		}
		return ret;
	}

	// -----------------------------
	// Non-local RM: Do the optimization

	// -------- Create new vector spaces without the nullspace
	giss::IndexTranslator trans_1_1p("trans_1_1p");
	trans_1_1p.init(n1(), used1);
	int n1p = trans_1_1p.nb();

	giss::IndexTranslator trans_3_3p("trans_3_3p");
	trans_3_3p.init(n3(), used3);
	int n3p = trans_3_3p.nb();

printf("n1p=%d, n3p=%d\n", n1p, n3p);

	// Allocate optimization problem
//	galahad::qpt_problem_c qpt(n1p, n3p, true, 1);	// Hessian_kind = 1
	galahad::qpt_problem_c qpt(n1p, n3p, true);	// m, n, eqp

	// ================== Constraints: RM f3 = f1
	// qpt.A = constraints matrix = RMp

	// LHS of constraints = RMp
	// Translate RM matrix while copying
	qpt.alloc_A(RM0->size());
	giss::ZD11SparseMatrix A_zd11(qpt.A, 0);
	for (auto ii = RM0->begin(); ii != RM0->end(); ++ii) {
		int i1 = ii.row();
		int i3 = ii.col();

		A_zd11.add(
			trans_1_1p.a2b(i1),
			trans_3_3p.a2b(i3),
			ii.val());
	}

	// Constraints: Ax + C = 0
	// Constraints RHS = f1p (rescaled)
	for (int i1p=0; i1p<n1p; ++i1p) {
		int i1 = trans_1_1p.b2a(i1p);
		qpt.C[i1p] = -f1(i1) * sum1_inv[i1];
	}

	// ================== Objective Function
	// Remember: qpt.Hessian_kind = -1;

	// H = 2 I
	qpt.alloc_H(n3p);
	giss::ZD11SparseMatrix H_zd11(qpt.H, 0);
	for (int i3p=0; i3p<n3p; ++i3p) {
		H_zd11.add(i3p, i3p, 2.0d);
	}
printf("CC1\n");
	// Fill in qpt.X0 = \Lambda f1
	qpt.f = 0;
	for (int i3p=0; i3p<n3p; ++i3p) {
		int i3 = trans_3_3p.b2a(i3p);

		int i1, k;
		hc_index->index_to_ik(i3, i1, k);
		double hval = f1(i1) * sum1_inv[i1];
//		qpt.X0[i3p] = hval;
		qpt.X[i3p] = hval;		// Starting value
		qpt.G[i3p] = -2.0d * hval;		// G = -2X0
		qpt.f += hval * hval;			// f = x0 . x0
	}

printf("CC2\n");

	// Solve it!
	double infinity = 1e20;
	eqp_solve_simple(qpt.this_f, infinity);

	// --------- Pick out the answer and convert back to standard vector space
	giss::CooVector<int, double> ret3;		// Function return value
	for (int i3p=0; i3p<n3p; ++i3p) {
		int i3 = trans_3_3p.b2a(i3p);
		ret3.add(i3, qpt.X[i3p]);
	}
	return ret3;
}
// --------------------------------------------------------------
/** TODO: This doesn't account for spherical earth */
std::unique_ptr<giss::VectorSparseMatrix> MatrixMaker::hp_to_atm()
{
//	int n1 = grid1->ndata();
printf("BEGIN hp_to_atm() %d %d\n", n1(), n3());
	std::unique_ptr<giss::VectorSparseMatrix> ret(
		new giss::VectorSparseMatrix(
		giss::SparseDescr(n1(), n3())));

	// Compute the hp->ice and ice->hc transformations for each ice sheet
	// and combine into one hp->hc matrix for all ice sheets.
	giss::SparseAccumulator<int,double> area1_m;
	for (auto sheet = sheets.begin(); sheet != sheets.end(); ++sheet) {
		auto hp2proj(sheet->hp_to_projatm(area1_m));
		if (correct_area1) hp2proj = multiply(
			*sheet->atm_proj_correct(ProjCorrect::PROJ_TO_NATIVE),
			*hp2proj);
		ret->append(*hp2proj);
	}

	giss::SparseAccumulator<int,double> area1_m_inv;
	divide_by(*ret, area1_m, area1_m_inv);
	ret->sum_duplicates();

printf("END hp_to_atm()\n");
	return ret;
}
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
// ==============================================================
// Write out the parts that this class computed --- so we can test/check them

boost::function<void ()> MatrixMaker::netcdf_define(NcFile &nc, std::string const &vname) const
{
	std::vector<boost::function<void ()>> fns;
	fns.reserve(sheets.size() + 1);

printf("MatrixMaker::netcdf_define(%s) (BEGIN)\n", vname.c_str());

	// ------ Attributes
	auto one_dim = giss::get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
	info_var->add_att("hc_index_type", _hc_index_type.str());

	// Names of the ice sheets
	std::string sheet_names = "";
	for (auto sheet = sheets.begin(); ; ) {
		sheet_names.append(sheet->name);
		++sheet;
		if (sheet == sheets.end()) break;
		sheet_names.append(",");
	}
	info_var->add_att("sheetnames", sheet_names.c_str());
#if 0
		info_var->add_att("grid1.name", gcm->grid1->name.c_str());
		info_var->add_att("grid2.name", grid2->name.c_str());
		info_var->add_att("exgrid.name", exgrid->name.c_str());
#endif

	// Define the variables
	fns.push_back(grid1->netcdf_define(nc, vname + ".grid1"));
	if (mask1.get())
		fns.push_back(giss::netcdf_define(nc, vname + ".mask1", *mask1));
	fns.push_back(giss::netcdf_define(nc, vname + ".hpdefs", hpdefs));
	for (auto sheet = sheets.begin(); sheet != sheets.end(); ++sheet) {
		fns.push_back(sheet->netcdf_define(nc, vname + "." + sheet->name));
	}


printf("MatrixMaker::netcdf_define(%s) (END)\n", vname.c_str());

	return boost::bind(&giss::netcdf_write_functions, fns);
}
// -------------------------------------------------------------
static std::vector<std::string> parse_comma_list(std::string list)
{
	std::stringstream ss(list);
	std::vector<std::string> result;

	while( ss.good() ) {
		std::string substr;
		getline( ss, substr, ',' );
		result.push_back( substr );
	}
	return result;
}

std::unique_ptr<IceSheet> read_icesheet(NcFile &nc, std::string const &vname)
{
	auto info_var = nc.get_var((vname + ".info").c_str());
	std::string stype(giss::get_att(info_var, "parameterization")->as_string(0));

	std::unique_ptr<IceSheet> sheet;
	if (stype == "L0") {
		sheet.reset(new IceSheet_L0);
	}
#if 0
	else if (stype == "L1") {
		sheet.reset(new IceSheet_L1);
	}
#endif

	sheet->read_from_netcdf(nc, vname);
	printf("read_icesheet(%s) END\n", vname.c_str());
	return sheet;

}


void MatrixMaker::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	clear();

	printf("MatrixMaker::read_from_netcdf(%s) 1\n", vname.c_str());
	grid1.reset(read_grid(nc, vname + ".grid1").release());
	if (giss::get_var_safe(nc, vname + ".mask1")) {
		mask1.reset(new blitz::Array<int,1>(
		giss::read_blitz<int,1>(nc, vname + ".mask1")));
printf("2 Set mask1 = %p\n", mask1.get());
	}
	hpdefs = giss::read_vector<double>(nc, vname + ".hpdefs");

	printf("MatrixMaker::read_from_netcdf(%s) 2\n", vname.c_str());

//	grid2.reset(read_grid(nc, "grid2").release());
//	exgrid.reset(read_grid(nc, "exgrid").release());

	// Read list of ice sheets
	NcVar *info_var = nc.get_var((vname + ".info").c_str());

	std::string shc_index_type(giss::get_att(info_var, "hc_index_type")->as_string(0));
	_hc_index_type = giss::parse_enum<HCIndex::Type>(shc_index_type.c_str());


	std::vector<std::string> sheet_names = parse_comma_list(std::string(
		giss::get_att(info_var, "sheetnames")->as_string(0)));

	for (auto sname = sheet_names.begin(); sname != sheet_names.end(); ++sname) {
		std::string var_name(vname + "." + *sname);
		printf("MatrixMaker::read_from_netcdf(%s) %s 3\n",
			vname.c_str(), var_name.c_str());
		add_ice_sheet(read_icesheet(nc, var_name));
	}

	// Remove grid cells that are not part of this domain.
	// TODO: This should be done while reading the cells in the first place.
	boost::function<bool (int)> include_cell1(domain->get_in_halo2());
	grid1->filter_cells(include_cell1);

	// Now remove cells from the exgrids and grid2s that interacted with grid1
	for (auto sheet=sheets.begin(); sheet != sheets.end(); ++sheet) {
		sheet->filter_cells1(include_cell1);
	}

}

std::unique_ptr<IceSheet> new_ice_sheet(Grid::Parameterization parameterization)
{
	switch(parameterization.index()) {
		case Grid::Parameterization::L0 : {
			IceSheet *ics = new IceSheet_L0;
			return std::unique_ptr<IceSheet>(ics);
//			return std::unique_ptr<IceSheet>(new IceSheet_L0);
		} break;
#if 0
		case Grid::Parameterization::L1 :
			return std::unique_ptr<IceSheet>(new IceSheet_L1);
		break;
#endif
		default :
			fprintf(stderr, "Unrecognized parameterization: %s\n", parameterization.str());
			throw std::exception();
	}
}


}
