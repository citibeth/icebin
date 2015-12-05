
#if 0
// Supports up to 64 ice sheets with 32-bit ints
// This can be expanded by moving from int to std::uint64 for index type in SparseVector and SpareMatrix
int const iI_mult = 33554432;

/** Gather index */
long i2_to_iI(int sheetno, int i2)
	{ return sheetno * iI_mult + i2; }

/** Scatter index */
void iI_to_i2(int iI, int &sheetno, int &i2) {
	shetno = iI / iI_mult;
	i2 = iI - sheetno * iI_mult;
}


/** Gathers vectors on each ice sheet into one vector for all ice sheets */
void MatrixMaker::gather_ice(
std::vector<std::pair<int, blitz::Array<double,1>>> const &srcs,
giss::MapSparseVector<int,double> &dest)
{
	dest.clear();
	for (auto ii=srcs.begin(); ii != srcs.end(); ++ii) {
		int sheet_id = srcs->first;
		//IceSheet *sheet = (*this)[srcs->first];
		blitz::Array<double,1> &src(srcs->second);

		size_t n2 = sheet->n2();
		for (int i2=0; i2<n2; ++i2) {
			if (sheet->masked(i2)) continue;

			long iI = sheet_id * iI_mult + i2;
			dest(iI) = i2;
		}
	}
}

#endif
