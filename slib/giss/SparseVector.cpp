void VectorSparseVector::sum_duplicates()
{
	sort_stable();

	// Scan through, overwriting our array
	// New array will never be bigger than original
	int j=-1;		// Last-written item in output array
	double last_index = -1;	// Last index we saw in input array
	for (int i=0; i<this->size(); ++i) {
		int index = (*this)[i].first;
		double val = (*this)[i].second;

		if (index == last_index) {
			(*this)[j].second += val;
		} else {
			++j;
			(*this)[j].second = val;
			last_index = index;
		}
	}
	int n = j+1;	// Size of output array

	this->resize(n);
	this->shrink_to_fit();
}
// ---------------------------------------------------
void VectorSparseVector::remove_duplicates()
{
	sort_stable();

	// Scan through, overwriting our array
	// New array will never be bigger than original
	int j=-1;		// Last-written item in output array
	int last_index = -1;	// Last index we saw in input array
	for (int i=0; i<this->size(); ++i) {
		int index = (*this)[i].first;
		double val = (*this)[i].second;

		if (index == last_index) {
			(*this)[j].second == val;
		} else {
			++j;
			(*this)[j].second = val;
			last_index = index;
		}
	}
	int n = j+1;	// Size of output array

	this->resize(n);
	this->shrink_to_fit();
}
