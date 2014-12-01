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

namespace glint2 {

static double const nan = std::numeric_limits<double>::quiet_NaN();

// REMEMBER: Decoding converts a set of (index, value) pairs into
// normal arrays (with NaN where no value was given.)
void IceModel_Decode::run_timestep(double time_s,
	blitz::Array<int,1> const &indices,
	std::vector<blitz::Array<double,1>> const &ivals2,
	std::vector<blitz::Array<double,1>> &ovals2)			// Output variables; should be allocated by caller
{
printf("BEGIN IceModel_Decode::run_timestep(time_s = %f) size=%ld\n", time_s, indices.size());
	std::vector<blitz::Array<double,1>> ivals2d;	/// Decoded fields

	// Naming convention on array variables:
	//     ivals2 = Vector of Values-arrays on grid2 (ice grid)
	//     ivals2d = Vector of DECODED values-arrays on grid2
	//     vals = Individual value array from ivals2
	//     valsd = Individual valu array from ivals2d
	// Loop through the fields we require
	int i=0;
	for (auto ii = ivals2.begin(); ii != ivals2.end(); ++ii, ++i) {
		blitz::Array<double,1> const &vals(*ii);

		// Decode the field!
		blitz::Array<double,1> valsd(ndata());
		valsd = nan;
		int n = indices.size();
		for (int i=0; i < n; ++i) {
			int ix = indices(i);
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

		// Store decoded field in our output
		ivals2d.push_back(valsd);
printf("Done decoding required field, %s\n", contract[IceModel::INPUT].name(i).c_str());
	}

	// Pass decoded fields on to subclass
	run_decoded(time_s, ivals2d, ovals2);
printf("END IceModel_Decode::run_timestep(%f)\n", time_s);
}


}
