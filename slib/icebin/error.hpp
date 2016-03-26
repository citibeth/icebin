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

#ifndef ICEBIN_ERROR_HPP
#define ICEBIN_ERROR_HPP

/** @defgroup ibmisc ibmisc.hpp
@brief Basic stuff common to all ibmisc */
namespace icebin {

// Use this instead.
// http://www.thecodingforums.com/threads/function-pointers-to-printf.317925/
/** @brief Printf-like signature of error handle functions to be used by SpSparse. */
typedef void (*error_ptr) (int retcode, char const *str, ...);

/** @brief Error handler used by IBMisc.  May be changed by user's
main program, to fit into some larger error handling system (eg:
Everytrace).

https://github.com/citibob/everytrace */
extern error_ptr icebin_error;

/** @brief Excpetion thrown by the default SpSparse error handler. */
class Exception : public std::exception
{
public:
    virtual ~Exception()
        {}

    virtual const char* what() const noexcept
        { return "icebin::Exception()"; }
};

}
/** @} */

#endif // Guard
