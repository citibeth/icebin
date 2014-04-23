#include <mpi.h>		// Must be first
#include <glint2/GCMParams.hpp>

namespace glint2 {

GCMParams::GCMParams() :
	gcm_rank(-1), gcm_root(-1) {}

void GCMParams::set_start_time(
	giss::time::tm const &_time_base,
	double _time_start_s)
{
	time_base = _time_base;
	time_start_s = _time_start_s;
}

GCMParams::GCMParams(
	MPI_Comm const _gcm_comm,
	int _gcm_root,
	boost::filesystem::path const &_config_dir)
: gcm_comm(_gcm_comm), gcm_root(_gcm_root), config_dir(_config_dir), time_start_s(0)
{
	MPI_Comm_rank(gcm_comm, &gcm_rank);

}

}

