

NcDim getOrAddDim(NcFile &nc, std::string const &dim_name, size_t dim_size = UNLIMITED)
{
	NcDim dim;
	try {
		// Look for the dim already existing
		dim = nc.getDim(dim_name);

		// Check that our request matches the existing dimension
		if (dim_size == UNLIMITED) {
			// Make sure existing dimension is unlimited
			if (!isUnlimited()) {
				std::string msg <<
					"Attempt in NcGroup::getOrAddDim to change size from " <<
					dim.getSize() << " to unlimited";
				throw NcBadDim(msg.c_str(), __FILE__, __LINE__);
			}
		} else {
			size_t existing_size = dim.getSize();
			if (existing_size != dim_size) {
				std::string msg <<
					"Attempt in NcGroup::getOrAddDim to change size from " <<
					existing_size << " to " << dim_size;
				throw NcBadDim(msg.c_str(), __FILE__, __LINE__);
			}
		}
	} except exceptions::NcBadDim {
		// Must create the new dim
		dim = nc.addDim(dim_name, dim_size);
	}

	return dim;
}
