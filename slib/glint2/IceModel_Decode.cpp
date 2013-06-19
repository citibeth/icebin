#include <glint2/IceModel_Decode.hpp>

namespace glint2 {

static double const nan = std::numeric_limits<double>::quiet_NaN();

void IceModel_Decode::run_timestep(long itime,
	blitz::Array<int,1> const &indices,
	std::map<IceField, blitz::Array<double,1>> const &vals2)
{
printf("BEGIN IceModel_Decode::run_timestep(%ld)\n", itime);
	std::map<IceField, blitz::Array<double,1>> vals2d;	/// Decoded fields

	// Loop through the fields we require
	std::set<IceField> fields;
printf("CC0\n");
	get_required_fields(fields);
printf("CC0\n");
	for (auto field = fields.begin(); field != fields.end(); ++field) {
printf("CC1a\n");
printf("Looking for required field %s\n", field->str());

		// Look up the field we require
		auto ii = vals2.find(*field);
		if (ii == vals2.end()) {
			fprintf(stderr, "Cannot find required ice field = %s\n", field->str());
			throw std::exception();
		}
		blitz::Array<double,1> vals(ii->second);

printf("CC1 ndata=%d\n", ndata);
		// Decode the field!
		blitz::Array<double,1> valsd(ndata);
		valsd = nan;
		int n = indices.size();
		for (int i=0; i < n; ++i) {
			int ix = indices(i);
			// Do our own bounds checking!
			if (ix < 0 || ix >= ndata) {
				fprintf(stderr, "IceModel: index %d out of range [0, %d)\n", ix, ndata);
				throw std::exception();
			}

#if 0
			// Sanity check for NaN coming through
			if (std::isnan(vals(i))) {
				fprintf(stderr, "IceModel::decode: vals[%d] (index=%d) is NaN!\n", i, ix);
				throw std::exception();
			}
#endif

			// Add this value to existing field
			double &oval = valsd(ix);
			if (std::isnan(oval)) oval = vals(i);
			else oval += vals(i);
		}
printf("CC2\n");

		// Store decoded field in our output
		vals2d.insert(std::make_pair(*field, valsd));
printf("CC3\n");
	}
printf("CC4\n");

	// Pass decoded fields on to subclass
	run_decoded(itime, vals2d);
printf("END IceModel_Decode::run_timestep(%ld)\n", itime);
}


}
