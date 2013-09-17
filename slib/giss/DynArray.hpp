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

#pragma once

namespace giss {

template<class T>
class DynArray {
public:
	size_t const size;
	size_t const ele_size;
	char * const buf;
	char * const buf_end;

	DynArray(size_t _ele_size, size_t _size) :
		ele_size(_ele_size),
		size(_size),
		buf((char *)malloc(size * ele_size)),
		buf_end(buf + size * ele_size) {}

	~DynArray() { free(buf); }

	DynArray(DynArray const &b) = delete;
	void operator=(DynArray const &b) = delete;
	void operator=(DynArray &&b) = delete;

	T *begin() { return reinterpret_cast<T *>(buf); }
	T *end() { return  reinterpret_cast<T *>(buf_end); }

	T &operator[](int ix)
		{ return *reinterpret_cast<T *>(buf + ele_size * ix); }

	T const &operator[](int ix) const
		{ return *reinterpret_cast<T const *>(buf + ele_size * ix); }

	void incr(T *&ptr) const {
		ptr = reinterpret_cast<T *>(
			reinterpret_cast<char *>(ptr) + ele_size);
	}

	void decr(T *&ptr) const {
		ptr = reinterpret_cast<T *>(
			reinterpret_cast<char *>(ptr) - ele_size);
	}

	/** Computes a-b */
	size_t diff(T *a, T*b) {
		size_t bdiff = reinterpret_cast<char *>(a) - reinterpret_cast<char *>(b);
		return bdiff / ele_size;
	}
};

}
