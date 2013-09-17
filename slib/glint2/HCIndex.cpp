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

#include <glint2/HCIndex.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/modele/ModelEDomain.hpp>

namespace glint2 {

std::unique_ptr<HCIndex> HCIndex::new_HCIndex(
	Type const type,
	MatrixMaker const &mm)
{
	switch(type.index()) {
		case HCIndex::Type::MODELE :
			return std::unique_ptr<HCIndex>(
				new glint2::modele::HCIndex_ModelE(mm.n1()));
	}
}

}
