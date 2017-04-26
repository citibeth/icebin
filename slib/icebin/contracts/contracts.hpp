/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <string>
#include <map>

namespace icebin {

class GCMCoupler;
class IceCoupler;

namespace contracts {

/** This field is returned at initialization time, before the first coupling. */
const unsigned INITIAL = 1;
const unsigned PRIVATE = 2;

extern std::string to_str(unsigned int flags);

// =====================================================
// Virtual Functions that dispatch on GCMCoupler and IceCoupler

void setup(GCMCoupler const *gcm_coupler, IceCoupler *ice_coupler);
// ======================================================



}}      // namespace
