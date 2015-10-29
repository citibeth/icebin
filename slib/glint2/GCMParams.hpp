#pragma once

#include <boost/filesystem.hpp>
#include <giss/time.hpp>

namespace glint2 {

/** Parameters passed from the GCM through to the ice model.
These parameters cannot be specific to either the ice model or the GCM. */
struct GCMParams {
	MPI_Comm gcm_comm;		// MPI communicator used by the GCM
	int gcm_rank;			// Rank of this process in gcm_comm
	int gcm_root;			// Rank of root process in gcm_comm
	boost::filesystem::path config_dir;	// Where to look for Ice Model configuration files
	boost::filesystem::path run_dir;		// Where to put output files
	giss::time::tm time_base;	// Corresponds to time_s == 0
	std::string time_units;		// CF-compliant string describing the time units
	double time_start_s;		// Start of simulation, as far as ice model is concerned (seconds since time_base).

	GCMParams();

	GCMParams(
		MPI_Comm const _gcm_comm,
		int _gcm_root,
		boost::filesystem::path const &_config_dir,
		boost::filesystem::path const &_run_dir);

	void set_start_time(
		giss::time::tm const &_time_base,
		double time_start_s);
};


}
