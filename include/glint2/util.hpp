#pragma once

namespace glint2 {

// See: http://stackoverflow.com/questions/37473/how-can-i-assert-without-using-abort
template <typename A>
inline void gassert(A assertion)
{
    if( !assertion ) throw std::exception();
}

}
