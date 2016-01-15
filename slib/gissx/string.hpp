#include <algorithm>
#include <string>

namespace giss {

inline void toupper(std::string &str) {
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

inline void tolower(std::string &str) {
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

inline bool ends_with(std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

}
