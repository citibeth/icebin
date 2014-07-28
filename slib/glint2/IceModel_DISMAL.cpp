/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013, 2014 by Robert Fischer
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
#include <boost/filesystem.hpp>
#include <glint2/IceModel_DISMAL.hpp>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <glint2/GCMParams.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {

void IceModel_DISMAL::IceModel_DISMAL::init(
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base)
{
	printf("BEGIN IceModel_DISMAL::init(%s)\n", vname_base.c_str());

	this->grid2 = grid2;

	// Transfer constants from GCM to PISM, and also set up coupling contracts.
	// This is the right place to do it, since the PISM systme is fully up and functional,
	// and all PISM config files have been read.
	// This call through the GCMCoupler will call back to setup_contracts_xxx().
	coupler->setup_contracts(*this);

	printf("END IceModel_DISMAL::int()\n");
}



/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_DISMAL::run_decoded(double time_s,
	std::vector<blitz::Array<double,1>> const &vals2)
{
}


}
