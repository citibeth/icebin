/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
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

#ifndef ICEBIN_SEARISE_INDEXING_HPP
#define ICEBIN_SEARISE_INDEXING_HPP

#include <icebin/Indexing.hpp>

namespace icebin {
namespace searise {


class Indexing : public icebin::Indexing {

public:

    const int nx, ny;

    Indexing(int _nx, int _ny) : nx(_nx), ny(_ny) {}


    long ij_to_index(int i, int j) const
        { return j*nx + i; }    // SEARISE ordering

    void index_to_ij(long index, int &i, int &j) const
    {
        j = index / nx;
        i = index - j*nx;
    }
};


}}

#endif  // Guard
