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

//#include <blitz/array.h>
//#include <giss/blitz.hpp>


#include <cstdio>
#include <vector>
#include <blitz/array.h>

namespace giss {

template<class T>
blitz::Array<T,1> vector_to_blitz(std::vector<T> &vec)
{
    blitz::TinyVector<int,1> shape(0);
    blitz::TinyVector<int,1> strides(0);

    shape[0] = vec.size();
    strides[0] = 1;     // Blitz++ strides in sizeof(T) units

    return blitz::Array<T,1>(&vec[0], shape, strides,
        blitz::neverDeleteData);
}

template<class T>
blitz::Array<T,1> const vector_to_blitz(std::vector<T> const &vec)
{
    blitz::TinyVector<int,1> shape(0);
    blitz::TinyVector<int,1> strides(0);

    shape[0] = vec.size();
    strides[0] = 1;     // Blitz++ strides in sizeof(T) units

    // const_cast because Blitz++ can't construct a const Array
    T *vecp = const_cast<T *>(&vec[0]);
    return blitz::Array<T,1>(vecp, shape, strides,
        blitz::neverDeleteData);
}

/** Makes sure a blitz::Array dimension.
Raises a Python exception if it does not. */
template<class T, int rank>
void check_dimensions(
std::string const &vname,
blitz::Array<T, rank> const &arr,
std::vector<int> const &dims)
{
    for (int i=0; i<rank; ++i) {
        if (dims[i] >= 0 && arr.extent(i) != dims[i]) {
            fprintf(stderr,
                "Error in %s: expected dimension #%d = %d (is %d instead)\n",
                vname.c_str(), i, dims[i], arr.extent(i));
            giss::exit(1);
        }
    }
}
// ------------------------------------------------------------
template<class T, int rank>
blitz::Array<T, rank> c_to_f(blitz::Array<T, rank> &arr);

template<class T, int rank>
blitz::Array<T, rank> c_to_f(blitz::Array<T, rank> &arr)
{
    // Initialize an 11-dim vector of transpositions
    // (because transpose() doesn't take a TinyVector)
    int const max_dims = 11;
    int rev[max_dims];
    for (int i=rank; i<max_dims; ++i) rev[i] = 0;

    // Reverse dimensions
    for (int i=0; i<rank; ++i) rev[i] = rank-1-i;
for (int i=0; i<11; ++i) printf("rev[%d] = %d\n", i, rev[i]);
    auto ret(arr.transpose(rev[0], rev[1], rev[2], rev[3], rev[4], rev[5], rev[6], rev[7], rev[8], rev[9], rev[10]));

    // Re-base to 1
    blitz::TinyVector<int, rank> base(1);
    ret.reindexSelf(base);

    return ret;
}
// ------------------------------------------------------------
template<class T, int rank>
blitz::Array<T, rank> f_to_c(blitz::Array<T, rank> &arr);

template<class T, int rank>
blitz::Array<T, rank> f_to_c(blitz::Array<T, rank> &arr)
{
    // Initialize an 11-dim vector of transpositions
    // (because transpose() doesn't take a TinyVector)
    int const max_dims = 11;
    int rev[max_dims];
    for (int i=rank; i<max_dims; ++i) rev[i] = 0;

    // Reverse dimensions
    for (int i=0; i<rank; ++i) rev[i] = rank-i;
    auto ret(arr.transpose(rev[0], rev[1], rev[2], rev[3], rev[4], rev[5], rev[6], rev[7], rev[8], rev[9], rev[10]));

    // Re-base to 0
    blitz::TinyVector<int, rank> base(0);
    ret.reindexSelf(base);

    return ret;
}

}

// ==============================================================

int main(int argc, char **argv)
{
    int nn[] = {2,3,4};

    blitz::Array<double,3> arr_c(nn[0], nn[1], nn[2]);
    arr_c = -1.0;
#if 1
    auto arr_f(giss::c_to_f(arr_c));
#else
    auto arr_f(arr_c.transpose(2,1,0,0,0,0,0,0,0,0,0));

    blitz::TinyVector<int, 3> base(1);
    arr_f.reindexSelf(base);
#endif

    for (int i=0; i<3; ++i) printf("stride %d: %d %d\n", i, arr_c.stride(i), arr_f.stride(i));

    for (int i=0; i<nn[0]; ++i) {
    for (int j=0; j<nn[1]; ++j) {
    for (int k=0; k<nn[2]; ++k) {
        printf("%p %p %f %f\n", &arr_c(i,j,k), &arr_f(k+1,j+1,i+1), arr_c(i,j,k), arr_f(k+1,j+1,i+1));
    }}}
}
