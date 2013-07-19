#include <giss/CooVector.hpp>
#include <giss/ncutil.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <glint2/HCIndex.hpp>
#include <giss/IndexTranslator.hpp>
#include <giss/IndexTranslator2.hpp>
#include <galahad/qpt_c.hpp>
#include <galahad/eqp_c.hpp>

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
/** @params f2 Some field on each ice grid (referenced by ID).  Do not have to be complete.
TODO: This only works on one ice sheet.  Will need to be extended
for multiple ice sheets. */
giss::CooVector<int, double>
MatrixMaker::ice_to_hp(
std::map<int, blitz::Array<double,1>> &f2s,
blitz::Array<double,1> &initial)
{
	// =============== Set up basic vector spaces for optimization problem
	std::set<int> used1, used3;
	std::set<std::pair<int,int>> used2;

	// Used in constraints
	std::unique_ptr<giss::VectorSparseMatrix> RM0(hp_to_atm());	// 3->1
	for (auto ii = RM0->begin(); ii != RM0->end(); ++ii) {
		used1.insert(ii.row());
		used3.insert(ii.col());
	}

	std::unique_ptr<giss::VectorSparseMatrix> RM(
		remove_small_constraints(*RM0, 2));
	RM0.reset();


	giss::SparseAccumulator<int,double> area1;
	giss::MapDict<int, giss::VectorSparseMatrix> Ss;
	giss::MapDict<int, giss::VectorSparseMatrix> XMs;
	std::map<int, size_t> size2;	// Size of each ice vector space
printf("f2s.size() = %d\n", f2s.size());
	for (auto f2i=f2s.begin(); f2i != f2s.end(); ++f2i) {
		IceSheet *sheet = (*this)[f2i->first];

		std::unique_ptr<giss::VectorSparseMatrix> S(
			sheet->ice_to_atm(area1));		// 2 -> 1
printf("S->size() = %d\n", S->size());
		for (auto ii = S->begin(); ii != S->end(); ++ii) {
			used1.insert(ii.row());
//printf("used.insert-a: %d %d\n", sheet->index, ii.col());
			used2.insert(std::make_pair(sheet->index, ii.col()));
		}

		std::unique_ptr<giss::VectorSparseMatrix> XM(
			sheet->hp_to_ice());				// 3 -> 2
		for (auto ii = XM->begin(); ii != XM->end(); ++ii) {
//printf("used.insert-b: %d %d\n", sheet->index, ii.row());
			used2.insert(std::make_pair(sheet->index, ii.row()));
			used3.insert(ii.col());
		}

		size2[sheet->index] = sheet->n2();

		// Store away for later reference
		Ss.insert(sheet->index, std::move(S));
		XMs.insert(sheet->index, std::move(XM));
	}

	giss::IndexTranslator trans_1_1p("trans_1_1p");
		trans_1_1p.init(n1(), used1);
printf("used2.size() = %d\n", used2.size());
	giss::IndexTranslator2 trans_2_2p("trans_2_2p");
		trans_2_2p.init(std::move(size2), used2);
	giss::IndexTranslator trans_3_3p("trans_3_3p");
		trans_3_3p.init(n3(), used3);

	int n1p = trans_1_1p.nb();
	int n2p = trans_2_2p.nb();
	int n3p = trans_3_3p.nb();

	// Translate to new matrices
	giss::VectorSparseMatrix RMp(giss::SparseDescr(n1p, n3p));
	giss::VectorSparseMatrix Sp(giss::SparseDescr(n1p, n2p));
	giss::VectorSparseMatrix XMp(giss::SparseDescr(n2p, n3p));

printf("n1p=%d, n2p=%d, n3p=%d\n", n1p, n2p, n3p);

printf("Translating RM\n");
	for (auto ii = RM->begin(); ii != RM->end(); ++ii) {
		RMp.add(
			trans_1_1p.a2b(ii.row()),
			trans_3_3p.a2b(ii.col()), ii.val());
	}


	for (auto f2i=f2s.begin(); f2i != f2s.end(); ++f2i) {
		int const index = f2i->first;
		IceSheet *sheet = (*this)[index];

		giss::VectorSparseMatrix *S(Ss[index]);
		giss::VectorSparseMatrix *XM(XMs[index]);

printf("Translating S: %d\n", index);
		for (auto ii = S->begin(); ii != S->end(); ++ii) {
			Sp.add(
				trans_1_1p.a2b(ii.row()),
				trans_2_2p.a2b(std::make_pair(index, ii.col())),
				ii.val());
		}

printf("Translating XM: %d\n", index);
		for (auto ii = XM->begin(); ii != XM->end(); ++ii) {
			XMp.add(
				trans_2_2p.a2b(std::make_pair(index, ii.row())),
				trans_3_3p.a2b(ii.col()),
				ii.val());
		}
	}

	// -------- Translate f2 -> f2p
	// Ignore elements NOT listed in the translation
	blitz::Array<double,1> f2p(n2p);
	f2p = 0;
	for (int i2p = 0; i2p < n2p; ++i2p) {
		std::pair<int,int> const &a(trans_2_2p.b2a(i2p));
		int index = a.first;
		int i2 = a.second;
		f2p(i2p) = f2s[index](i2);
	}

	// ----------- Translate area1 -> area1p
	blitz::Array<double,1> area1p_inv(n1p);
	area1p_inv = 0;
	for (auto ii = area1.begin(); ii != area1.end(); ++ii) {
		int i1 = ii->first;
		int i1p = trans_1_1p.a2b(i1);
		area1p_inv(i1p) += ii->second;
	}
	for (int i1p=0; i1p<n1p; ++i1p) {
		if (area1p_inv(i1p) != 0) area1p_inv(i1p) = 1.0d / area1p_inv(i1p);
	}

	// ---------- Divide Sp by area1p to complete the regridding matrix
	for (auto ii = Sp.begin(); ii != Sp.end(); ++ii) {
		int i1p = ii.row();
		ii.val() *= area1p_inv(i1p);
	}

	// ========================================================
	// ========================================================

 	// ---------- Allocate the QPT problem
	// m = # constraints = n1p (size of atmosphere grid)
	// n = # variabeles = n3p
	galahad::qpt_problem_c qpt(n1p, n3p, true);

	// ================ Objective Function
	// 1/2 (A F_E - F_I)^2    where A = XM = (Ice->Exch)(Elev->Ice)
	// qpt%H = A^T A,    qpt%G = f_I \cdot A,        qpt%f = f_I \cdot f_I

	// -------- H = 2 * XMp^T XMp
	giss::VectorSparseMatrix XMp_T(giss::SparseDescr(XMp.ncol, XMp.nrow));
	transpose(XMp, XMp_T);
printf("XMp: (%d x %d) = %d elements\n", XMp.nrow, XMp.ncol, XMp.size());
printf("XMp_T: (%d x %d) = %d elements\n", XMp_T.nrow, XMp_T.ncol, XMp_T.size());
	std::unique_ptr<giss::VectorSparseMatrix> H(multiply(XMp_T, XMp));	// n3xn3
printf("H: (%d x %d) = %d elements\n", H->nrow, H->ncol, H->size());
	qpt.alloc_H(H->size());
	giss::ZD11SparseMatrix H_zd11(qpt.H, 0);
	for (auto ii = H->begin(); ii != H->end(); ++ii)
		H_zd11.add(ii.row(), ii.col(), 2.0d * ii.val());

	// -------- Linear term of obj function
	// G = -2*f2p \cdot XMp
	for (int i=0; i < qpt.n; ++i) qpt.G[i] = 0;
	for (auto ii = XMp.begin(); ii != XMp.end(); ++ii) {
		qpt.G[ii.col()] -= 2.0d * f2p(ii.row()) * ii.val();
	}

	// --------- Constant term of objective function
	// f = f2p \cdot f2p
	qpt.f = 0;
	for (int i2p=0; i2p<n2p; ++i2p) {
		qpt.f += f2p(i2p) * f2p(i2p);
	}

	// De-allocate...
	H.reset();
	XMp.clear();
	XMp_T.clear();

	// ============================ Constraints
	// RM x = Sp f2p

	// qpt.A = constraints matrix = RMp
	qpt.alloc_A(RMp.size());
	giss::ZD11SparseMatrix A_zd11(qpt.A, 0);
	copy(RMp, A_zd11);

	// qpt.C = equality constraints RHS = Sp * f2p
	for (int i=0; i<n1p; ++i) qpt.C[i] = 0;
	for (auto ii = Sp.begin(); ii != Sp.end(); ++ii) {
		int i1p = ii.row();		// Atm
		int i2p = ii.col();		// Ice
		qpt.C[i1p] += f2p(i2p) * ii.val();
	}

	// De-allocate
	RMp.clear();

	// =========================== Initial guess at solution
	for (int i3p=0; i3p<n3p; ++i3p) {
		int i3 = trans_3_3p.b2a(i3p);
		qpt.X[i3p] = initial(i3);
//		qpt.X[i3p] = 0;	// we have no idea
	}

	// =========================== Solve the Problem!
	double infinity = 1e20;
	eqp_solve_simple(qpt.this_f, infinity);



	// ========================================================
	// ========================================================

	// --------- Pick out the answer and convert back to standard vector space
	giss::CooVector<int, double> ret;
	for (int i3p=0; i3p<n3p; ++i3p) {
		int i3 = trans_3_3p.b2a(i3p);
		ret.add(i3, qpt.X[i3p]);
	}

	return ret;
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
		ret->append(*sheet->hp_to_atm(area1_m));
	}

	giss::SparseAccumulator<int,double> area1_m_inv;
	divide_by(*ret, area1_m, area1_m_inv);
printf("After divide_by: %ld %d\n", area1_m.size(), area1_m_inv.size());
	ret->sum_duplicates();

#if 0
printf("Writing hp2atm ret = %p\n", ret.get());
NcFile nc("hp2atm.nc", NcFile::Replace);
ret->netcdf_define(nc, "hp2atm")();
nc.close();
printf("Done Writing hp2hc ret = %p\n", ret.get());
#endif

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
		fns.push_back(giss::netcdf_define(nc, vname + "mask1", *mask1));
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
	}
	hpdefs = giss::read_vector<double>(nc, vname + ".hpdefs");

	printf("MatrixMaker::read_from_netcdf(%s) 2\n", vname.c_str());

//	grid2.reset(read_grid(nc, "grid2").release());
//	exgrid.reset(read_grid(nc, "exgrid").release());

	// Read list of ice sheets
	NcVar *info_var = nc.get_var((vname + ".info").c_str());
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
