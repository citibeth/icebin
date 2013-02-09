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

