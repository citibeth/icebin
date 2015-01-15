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

#include <mpi.h>		// Must be first
#include <glint2/IceModel_Decode.hpp>
#include <giss/sort.hpp>

using namespace giss;

namespace glint2 {

static double const nan = std::numeric_limits<double>::quiet_NaN();

// REMEMBER: Decoding converts a set of (index, value) pairs into
// normal arrays (with NaN where no value was given.)
void IceModel_Decode::run_timestep(double time_s,
	blitz::Array<int,1> const &indices,
	std::vector<blitz::Array<double,1>> const &ivals2)
{
printf("BEGIN IceModel_Decode::run_timestep(time_s = %f) size=%ld\n", time_s, indices.size());

	blitz::Array<int,1> nindices;
	std::vector<blitz::Array<double,1>> nivals2;

	blitz::Array<int,1> const *xindices;
	std::vector<blitz::Array<double,1>> const *xivals2;

	// Test out freestanding sorting/consolidation code
	// This section of code is NOT necessary.
	if (false) {
		std::vector<int> perm = sorted_perm(indices);
		nindices.reference(blitz::Array<int,1>(indices.size()));
		int nconsolidated = consolidate_by_perm(indices, perm, indices, nindices, DuplicatePolicy::REPLACE);
printf("IceModel_Decode: consolidated from %d down to %d\n", indices.size(), nconsolidated);
		nindices.reference(blitz::Array<int,1>(nconsolidated));
		consolidate_by_perm(indices, perm, indices, nindices, DuplicatePolicy::REPLACE);

		for (unsigned int i=0; i<ivals2.size(); ++i) {
			blitz::Array<double,1> nvals(nconsolidated);
			consolidate_by_perm(indices, perm, ivals2[i], nvals, DuplicatePolicy::ADD);
			nivals2.push_back(nvals);
		}

		xindices = &nindices;
		xivals2 = &nivals2;
	} else {
		xindices = &indices;
		xivals2 = &ivals2;
	}

	std::vector<blitz::Array<double,1>> ivals2d;	/// Decoded fields

	// Naming convention on array variables:
	//     ivals2 = Vector of Values-arrays on grid2 (ice grid)
	//     ivals2d = Vector of DECODED values-arrays on grid2
	//     vals = Individual value array from ivals2
	//     valsd = Individual valu array from ivals2d
	// Loop through the fields we require
	giss::CouplingContract const &icontract(contract[IceModel::INPUT]);
	for (int i=0; i<icontract.size_nounit(); ++i) {

		blitz::Array<double,1> const &vals((*xivals2)[i]);

		// Decode the field!
		blitz::Array<double,1> valsd(ndata());
		valsd = nan;
		int n = xindices->size();
		for (int i=0; i < n; ++i) {
			int ix = (*xindices)(i);
			// Do our own bounds checking!
			if (ix < 0 || ix >= ndata()) {
				fprintf(stderr, "IceModel: index %d out of range [0, %d)\n", ix, ndata());
				throw std::exception();
			}

			// Add this value to existing field
			double &oval = valsd(ix);
			if (std::isnan(oval)) oval = vals(i);
			else oval += vals(i);
		}

		// Convert any remaining nans to default value,
		// so we have a valid number everywhere.
		double default_value = icontract.field(i).default_value;
		for (int j=0; j<ndata(); ++j) {
			double &val(valsd(j));
			if (std::isnan(val)) val = default_value;
		}

		// Store decoded field in our output
		ivals2d.push_back(valsd);
printf("Done decoding required field, %s\n", icontract.name(i).c_str());
	}

	// Pass decoded fields on to subclass
	run_decoded(time_s, ivals2d);
printf("END IceModel_Decode::run_timestep(%f)\n", time_s);
}


}
