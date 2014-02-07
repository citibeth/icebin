#include <mpi.h>		// Must be first
#include <glint2/IceModel.hpp>

namespace glint2 {


IceModel::GCMParams::GCMParams() :
	gcm_rank(-1), gcm_root(-1) {}

IceModel::GCMParams::GCMParams(
	MPI_Comm const _gcm_comm,
	int _gcm_root,
	boost::filesystem::path const &_config_dir,
	giss::time::tm const &_time_base)
: gcm_comm(_gcm_comm), gcm_root(_gcm_root), config_dir(_config_dir), time_base(_time_base)
{
	MPI_Comm_rank(gcm_comm, &gcm_rank);

}


}
