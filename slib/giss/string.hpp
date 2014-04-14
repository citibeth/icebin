#include <algorithm>
#include <string>

namespace giss {

void toupper(std::string &str) {
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

void tolower(std::string &str) {
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

}
