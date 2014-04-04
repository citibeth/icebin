#include <giss/VarTransformer.hpp>

namespace giss {



void VarTransformer::set(std::string output, std::string input, std::string scalar, double val)
{
	int ioutput = name_to_ix(OUTPUTS, output);
	if (ioutput < 0) return;		// No error if our output variable is not in this tensor.
	int iinput = name_to_ix(INPUTS, input);
	int iscalar = name_to_ix(SCALARS, scalar);

	tensor(ioutput, iinput, iscalar) = val;
}


/** Instantiates the scalars with specific values, and returns a 2nd-order
matrix derived from the 3d-order tensor, in CSR format. */
CSRAndUnits VarTransformer::apply_scalars(
		std::vector<std::pair<std::string, double>> const &nvpairs)
{
	// Convert name/value pairs to a regular vector
	blitz::Array<double,1> scalars(dimsize(SCALARS));
	scalars = 0;
	for (auto ii = nvpairs.begin(); ii != nvpairs.end(); ++ii) {
		std::string const &name = ii->first;
		double const val = ii->second;
		scalars[name_to_ix(SCALARS, name)] = val;
	}

	// Take inner product of tensor with our scalars.
	CSRAndUnits ret(dimsize(OUTPUTS));
	for (int i=0; i < dimsize(OUTPUTS); ++i) {
		int unit_ix = dimsize(INPUTS)-1;
		for (int j=0; j < dimsize(INPUTS); ++i) {
			double coeff = 0;
			for (int k=0; k < dimsize(SCALARS); ++i) {
				coeff += tensor(i,j,k) * scalars(k);
			}

			// Output format: sparse matrix plus dense unit column
			if (j == unit_ix) {
				ret.units[j] = coeff;
			} else {
				if (coeff != 0) ret.mat.add(i, j, coeff);
			}
		}
	}

	return ret;
}

std::ostream &VarTransformer::operator<<(std::ostream &out)
{
	int unit_k = dimsize(SCALARS) - 1;
	int unit_j = dimsize(INPUTS) - 1;
	for (int i=0; i<dimsize(OUTPUTS); ++i) {
		out << (*ele_names[OUTPUTS])[i] << " = ";

		// Count number of INPUTs used for this OUTPUT
		int nj = 0;
		std::vector<int> nk(dimsize(INPUTS), 0);
		for (int j=0; j < dimsize(INPUTS); ++j) {
			for (int k=0; k < dimsize(SCALARS); ++k) {
				if (tensor(i,j,k) != 0) ++nk[j];
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
		for (int j=0; j < dimsize(INPUTS); ++j) {
			int nkj = nk[j];
			if (nkj == 0) continue;

			if (nkj > 1) out << "(";
			int kk = 0;
			for (int k=0; k < dimsize(SCALARS); ++k) {
				double val = tensor(i,j,k);
				if (val == 0.0) continue;
				out << val;
				if (k != unit_k) out << " " << (*ele_names[SCALARS])[k];

				if (kk != nkj-1) out << " + ";

				// Increment count of SEEN k values
				++kk;
			}
			if (nkj > 1) out << ")";
			if (j != unit_j) out << " " << (*ele_names[INPUTS])[j];

			if (jj != nj-1) out << " + ";

			out << (*ele_names[INPUTS])[j];

			// Increment count of SEEN j values
			++jj;
		}
	}
	return out;
}

}	// namespace giss