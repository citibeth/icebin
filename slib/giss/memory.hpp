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

// See: http://stackoverflow.com/questions/11002641/dynamic-casting-for-unique-ptr

#include <memory>
#include <typeinfo>

namespace giss {

template <class T_DEST, class T_SRC>
inline std::unique_ptr<T_DEST> unique_cast(std::unique_ptr<T_SRC> &&src)
{
	if (!src) return std::unique_ptr<T_DEST>();

	// Throws a std::bad_cast() if this doesn't work out
	T_DEST *dest_ptr = &dynamic_cast<T_DEST &>(*src.get());

	src.release();
	return std::unique_ptr<T_DEST>(dest_ptr);
}

template <class T_DEST, class T_SRC>
inline std::shared_ptr<T_DEST> shared_cast(std::unique_ptr<T_SRC> &&src)
{
	return std::shared_ptr<T_DEST>(unique_cast<T_DEST, T_SRC>(std::move(src)));
}

}	// namespace giss








#if 0
template <typename T_SRC, typename T_DEST, typename T_DELETER>
bool dynamic_pointer_move(std::unique_ptr<T_DEST, T_DELETER> & dest,
                          std::unique_ptr<T_SRC, T_DELETER> & src)
{
    if (!src) {
        dest.reset();
        return true;
    }

    T_DEST * dest_ptr = dynamic_cast<T_DEST *>(src.get());
    if (!dest_ptr)
        return false;

    std::unique_ptr<T_DEST, T_DELETER> dest_temp(dest_ptr, src.get_deleter());
    src.release();
    dest.swap(dest_temp);
    return true;
}

template <typename T_SRC, typename T_DEST>
bool dynamic_pointer_move(std::unique_ptr<T_DEST> & dest,
                          std::unique_ptr<T_SRC> & src)
{
    if (!src) {
        dest.reset();
        return true;
    }

    T_DEST * dest_ptr = dynamic_cast<T_DEST *>(src.get());
    if (!dest_ptr)
        return false;

    src.release();
    dest.reset(dest_ptr);
    return true;
}
#endif

