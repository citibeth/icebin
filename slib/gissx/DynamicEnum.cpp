#if 0
#include <giss/DynamicEnum.hpp>

namespace giss {

void DynamicEnum::set_names(std::vector &&names)
{
	_ix_to_name = std::move(names);
	_name_to_ix.clear();
	int i = 0;
	for (auto ii = _ix_to_name.begin(); ii != _ix_to_name.end(); ++ii, ++i) {
		_name_to_ix[dim][i] = *ii;
	}
}


}
#endif
