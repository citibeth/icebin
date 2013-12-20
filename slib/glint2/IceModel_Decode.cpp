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

void IceModel_Decode::run_timestep(double time_s,
	blitz::Array<int,1> const &indices,
	std::map<IceField, blitz::Array<double,1>> const &vals2)
{
printf("BEGIN IceModel_Decode::run_timestep(%f)\n", time_s);
	std::map<IceField, blitz::Array<double,1>> vals2d;	/// Decoded fields

	// Loop through the fields we require
	std::set<IceField> fields;
	get_required_fields(fields);
	for (auto field = fields.begin(); field != fields.end(); ++field) {
printf("Looking for required field %s\n", field->str());

		// Look up the field we require
		auto ii = vals2.find(*field);
		if (ii == vals2.end()) {
			fprintf(stderr, "Cannot find required ice field = %s\n", field->str());
			throw std::exception();
		}
		blitz::Array<double,1> vals(ii->second);

		// Decode the field!
		blitz::Array<double,1> valsd(ndata);
		valsd = nan;
		int n = indices.size();
		for (int i=0; i < n; ++i) {
			int ix = indices(i);
			// Do our own bounds checking!
			if (ix < 0 || ix >= ndata) {
				fprintf(stderr, "IceModel: index %d out of range [0, %d)\n", ix, ndata);
				throw std::exception();
			}

#if 0
			// Sanity check for NaN coming through
			if (std::isnan(vals(i))) {
				fprintf(stderr, "IceModel::decode: vals[%d] (index=%d) is NaN!\n", i, ix);
				throw std::exception();
			}
#endif

			// Add this value to existing field
			double &oval = valsd(ix);
			if (std::isnan(oval)) oval = vals(i);
			else oval += vals(i);
		}

		// Store decoded field in our output
		vals2d.insert(std::make_pair(*field, valsd));
	}

	// Pass decoded fields on to subclass
	run_decoded(time_s, vals2d);
printf("END IceModel_Decode::run_timestep(%ld)\n", time_s);
}


}
