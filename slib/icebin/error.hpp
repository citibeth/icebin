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
