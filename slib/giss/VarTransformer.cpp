#include <giss/VarTransformer.hpp>

namespace giss {

void VarTransformer::allocate()
{
	blitz::TinyVector<int, NDIM> extent;
	for (int i=0; i<NDIM; ++i) extent[i] = _ele_names[i]->size_withunit();
	_tensor.reference(blitz::Array<double, NDIM>(extent));
	_tensor = 0;
}


bool VarTransformer::set(std::string output, std::string input, std::string scalar, double val)
{
	int ioutput = dimension(OUTPUTS)[output];
	int iinput = dimension(INPUTS)[input];
	int iscalar = dimension(SCALARS)[scalar];

	bool ret = true;
	if (ioutput < 0) {
		fprintf(stderr, "ERROR: VarTransformer::set(): output variable %s not defined.\n", output.c_str());
		ret = false;
	}
	if (iinput < 0) {
		fprintf(stderr, "ERROR: VarTransformer::set(): input variable %s not defined.\n", input.c_str());
		ret = false;
	}

	if (iscalar < 0) {
		fprintf(stderr, "ERROR: VarTransformer::set(): scalar variable %s not defined.\n", output.c_str());
		ret = false;
	}

	if (!ret) return ret;
	
	_tensor(ioutput, iinput, iscalar) = val;
	return ret;
}


/** Instantiates the scalars with specific values, and returns a 2nd-order
matrix derived from the 3d-order tensor, in CSR format. */
CSRAndUnits VarTransformer::apply_scalars(
		std::vector<std::pair<std::string, double>> const &nvpairs)
{
	printf("BEGIN VarTransformer::apply_scalars()\n");
	for (auto ele : nvpairs) {
		printf("    %s = %g\n", ele.first.c_str(), ele.second);
	}


	int n_outputs_nu = dimension(OUTPUTS).size_nounit();		// # OUTPUTS no unit
	int n_inputs_wu = dimension(INPUTS).size_withunit();
	int n_scalars_wu = dimension(SCALARS).size_withunit();	// # SCALARS w/unit

	int unit_inputs = dimension(INPUTS).unit_ix();

	// Convert name/value pairs to a regular vector
	blitz::Array<double,1> scalars(n_scalars_wu);
	scalars = 0;
	for (auto ii = nvpairs.begin(); ii != nvpairs.end(); ++ii) {
		std::string const &nv_name = ii->first;
		double const val = ii->second;
		scalars(dimension(SCALARS)[nv_name]) = val;
	}

//std::cout << "Input vector = " << scalars << std::endl;

	// Take inner product of tensor with our scalars.
	CSRAndUnits ret(n_outputs_nu);
	for (int i=0; i < n_outputs_nu; ++i) {
		for (int j=0; j < n_inputs_wu; ++j) {
			double coeff = 0;
			for (int k=0; k < n_scalars_wu; ++k) {
				coeff += _tensor(i,j,k) * scalars(k);
			}

			// Output format: sparse matrix plus dense unit column
			if (j == unit_inputs) {
				ret.units[i] = coeff;
			} else {
				if (coeff != 0) ret.mat.add(i, j, coeff);
			}
		}
	}

	std::cout << "apply_scalars() returning " << ret;

	printf("END VarTransformer::apply_scalars()\n");
	return ret;
}

std::ostream &operator<<(std::ostream &out, CSRMatrix const &mat)
{
	out << "CSRMatrix :" << std::endl;
	for (unsigned int i=0; i < mat.size(); ++i) {
		out << "    " << i << ":";
		for (auto ele : mat[i]) {
			out << " (" << ele.first << ", " << ele.second << ")";
		}
		out << std::endl;
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, CSRAndUnits const &matu)
{
	CSRMatrix const &mat(matu.mat);

	out << "CSRMatrix :" << std::endl;
	for (unsigned int i=0; i < mat.size(); ++i) {
		out << "    " << i << ": " << matu.units[i] << " +";
		for (auto ele : mat[i]) {
			out << " (" << ele.first << ", " << ele.second << ")";
		}
		out << std::endl;
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, VarTransformer const &vt)
{
	int n_outputs_nu = vt.dimension(VarTransformer::OUTPUTS).size_nounit();		// # OUTPUTS no unit
	int n_inputs_wu = vt.dimension(VarTransformer::INPUTS).size_withunit();
	int n_scalars_wu = vt.dimension(VarTransformer::SCALARS).size_withunit();	// # SCALARS w/unit

	int unit_outputs = vt.dimension(VarTransformer::OUTPUTS).unit_ix();
	int unit_inputs = vt.dimension(VarTransformer::INPUTS).unit_ix();
	int unit_scalars = vt.dimension(VarTransformer::SCALARS).unit_ix();


	for (int i=0; i<n_outputs_nu; ++i) {
		out << "    " << vt.dimension(VarTransformer::OUTPUTS)[i] << " = ";

		// Count number of INPUTs used for this OUTPUT
		int nj = 0;
		std::vector<int> nk(n_inputs_wu, 0);
		for (int j=0; j < n_inputs_wu; ++j) {
			for (int k=0; k < n_scalars_wu; ++k) {
				if (vt._tensor(i,j,k) != 0) ++nk[j];
			}
			if (nk[j] > 0) ++nj;
		}

		// No RHS for this OUTPUT, quit
		if (nj == 0) {
			out << "0" << std::endl;
			continue;
		}

		// We DO have something on the RHS
		int jj = 0;
		for (int j=0; j < n_inputs_wu; ++j) {
			int nkj = nk[j];
			if (nkj == 0) continue;

			if (nkj > 1) out << "(";
			int kk = 0;
			for (int k=0; k < n_scalars_wu; ++k) {
				double val = vt._tensor(i,j,k);
				if (val == 0.0) continue;
				if (val != 1.0) out << val;
				if (k != unit_scalars) out << " " << vt.dimension(VarTransformer::SCALARS)[k];

				if (kk != nkj-1) out << " + ";

				// Increment count of SEEN k values
				++kk;
			}
			if (nkj > 1) out << ")";
			if (j != unit_inputs) out << " " << vt.dimension(VarTransformer::INPUTS)[j];

			if (jj != nj-1) out << " + ";

//			out << vt.dimension(VarTransformer::INPUTS)[j];

			// Increment count of SEEN j values
			++jj;
		}
		out << std::endl;
	}
	return out;
}

}	// namespace giss
