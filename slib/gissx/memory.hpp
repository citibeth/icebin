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

// http://ficksworkshop.com/blog/14-coding/86-how-to-static-cast-std-unique-ptr	
template<typename D, typename B>
std::unique_ptr<D> static_cast_unique_ptr(std::unique_ptr<B>& base)
{
    return std::unique_ptr<D>(static_cast<D*>(base.release()));
}
  
template<typename D, typename B>
std::unique_ptr<D> static_cast_unique_ptr(std::unique_ptr<B>&& base)
{
    return std::unique_ptr<D>(static_cast<D*>(base.release()));
}


template<typename D, typename B>
std::unique_ptr<D> dynamic_cast_unique_ptr(std::unique_ptr<B>& base)
{
    return std::unique_ptr<D>(dynamic_cast<D*>(base.release()));
}
  
template<typename D, typename B>
std::unique_ptr<D> dynamic_cast_unique_ptr(std::unique_ptr<B>&& base)
{
    return std::unique_ptr<D>(dynamic_cast<D*>(base.release()));
}



template <class T_DEST, class T_SRC>
inline std::shared_ptr<T_DEST> dynamic_shared_cast(std::unique_ptr<T_SRC> &&src)
{
	return std::shared_ptr<T_DEST>(dynamic_cast_unique_ptr<T_DEST, T_SRC>(std::move(src)));
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

