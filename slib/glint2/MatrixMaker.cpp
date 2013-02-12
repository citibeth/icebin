#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <glint2/eigen.hpp>
#include <giss/ncutil.hpp>

namespace glint2 {

void MatrixMaker::clear()
{
	sheets.clear();
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

	long nhc = hpdefs.size();
	if (hcmax.extent(0) != nhc-1) {
		fprintf(stderr, "hcmax for %s has wrong size: %d (vs %d expected)\n",
			mask1->extent(0), n1);
		throw std::exception();
	}

	// ------------- Realize the ice sheets
	for (auto ii=sheets.begin(); ii != sheets.end(); ++ii)
		(*ii)->realize();
}

int MatrixMaker::add_ice_sheet(std::unique_ptr<IceSheet> &&sheet) {
	sheet->gcm = this;
	int n = sheets.size();
	sheets.push_back(std::move(sheet));
	return n;
}


/** NOTE: Does not necessarily assume that ice sheets do not overlap on the same GCM grid cell */
void MatrixMaker::compute_fhc(
	blitz::Array<double,2> *fhc1h,	// OUT
	blitz::Array<double,1> *fgice1)	// OUT: Portion of gridcell covered in ground ice (from landmask)
{
	// Zero it all out
	int n1 = grid1->ndata();
	if (fhc1h) {
		for (int ihc=0; ihc<nhc(); ++ihc) {
			for (int i1=0; i1<n1; ++i1) {
				(*fhc1h)(ihc,i1) = 0.;
			}
		}
	}

	if (fgice1) {
		for (int i1=0; i1<n1; ++i1) (*fgice1)(i1) = 0.;
	}

	// Add in for each ice sheet
	for (auto sheet = sheets.begin(); sheet != sheets.end(); ++sheet) {
		(*sheet)->compute_fhc(fhc1h, fgice1);
	}
}

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
	for (auto sheetp = sheets.begin(); ; ) {
		IceSheet &sheet = **sheetp;
		sheet_names.append(sheet.name);
		++sheetp;
		if (sheetp == sheets.end()) break;
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
printf("***************** 1 hcmax.extent(0) = %d\n", hcmax.extent(0));
	fns.push_back(giss::netcdf_define(nc, vname + ".hcmax", hcmax));
	for (auto sheetp = sheets.begin(); sheetp != sheets.end(); ++sheetp) {
		IceSheet &sheet = **sheetp;
		fns.push_back(sheet.netcdf_define(nc, vname + "." + sheet.name));
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

std::unique_ptr<IceSheet> read_ice_sheet(NcFile &nc, std::string const &vname)
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
	printf("read_ice_sheet(%s) END\n", vname.c_str());
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
	hcmax.reference(giss::read_blitz<double,1>(nc, vname + ".hcmax"));

	printf("MatrixMaker::read_from_netcdf(%s) 2\n", vname.c_str());

//	grid2.reset(read_grid(nc, "grid2").release());
//	exgrid.reset(read_grid(nc, "exgrid").release());

	// Read list of ice sheets
	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	std::vector<std::string> sheet_names(parse_comma_list(std::string(
		giss::get_att(info_var, "sheetnames")->as_string(0))));

	for (auto sname = sheet_names.begin(); sname != sheet_names.end(); ++sname) {
		std::string sheet_name(vname + "." + *sname);
		printf("MatrixMaker::read_from_netcdf(%s) %s 3\n", vname.c_str(), sheet_name.c_str());
		sheets.push_back(read_ice_sheet(nc, sheet_name));
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
