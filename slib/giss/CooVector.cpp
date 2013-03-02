#include <giss/CooVector.hpp>

namespace giss {

struct CmpIndex1 {
	int *index1;
public:
	CmpIndex1(int *_index1) {
		index1 = _index1;
	}
	bool operator()(int i, int j)
	{
		return (index1[i] < index1[j]);
	}
};

void CooVector::sort()
{
	// Decide on how we'll sort
	CmpIndex1 cmp(&indices[0]);

	// Generate a permuatation
	int n = size();
	std::vector<int> perm; perm.reserve(n);
	for (int i=0; i<n; ++i) perm.push_back(i);
	std::sort(perm.begin(), perm.end(), cmp);

	// Apply permutation to vals
	std::vector<double> dtmp; dtmp.reserve(n);
	for (int i=0; i<n; ++i) dtmp.push_back(vals[perm[i]]);
	vals = std::move(dtmp);

	// Apply permutation to indices
	std::vector<int> itmp; itmp.reserve(n);
	for (int i=0; i<n; ++i) itmp.push_back(indices[perm[i]]);
	indices = std::move(itmp);
}

void CooVector::sum_duplicates()
{
	// Decide on how we'll sort
	CmpIndex1 cmp(&indices[0]);

	// Generate a sorted permuatation
	int n = size();
	std::vector<int> perm; perm.reserve(n);
	for (int i=0; i<n; ++i) perm.push_back(i);
	std::sort(perm.begin(), perm.end(), cmp);

	// Output arrays
	std::vector<int> nindices;
	std::vector<double> nvals;

	// Identify duplicates
	nindices.push_back(indices[perm[0]]);
	nvals.push_back(vals[perm[0]]);
	for (int i=1; i<indices.size(); ++i) {
		int row = indices[perm[i]];
//printf("remove_dup: %d %d\n", row, col);
		if ((row == nindices.back())) {
			nvals.back() += vals[perm[i]];
		} else {
			nindices.push_back(indices[perm[i]]);
			nvals.push_back(vals[perm[i]]);
		}
	}

	indices = std::move(nindices);
	vals = std::move(nvals);
}



}
