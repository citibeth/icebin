#include <glint2/HCIndex.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/modele/ModelEDomain.hpp>

namespace glint2 {

std::unique_ptr<HCIndex> HCIndex::new_HCIndex(
	Type const type,
	MatrixMaker const &mm)
{
	switch(type.index()) {
		case HCIndex::Type::MODELE :
			return std::unique_ptr<HCIndex>(
				new glint2::modele::HCIndex_ModelE(mm.n1()));
	}
}

}
