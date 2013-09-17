/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <netcdfcpp.h>
#include <giss/SparseMatrix.hpp>

namespace giss {


// -------------------------------------------------------
struct CmpIndex2 {
	int *index1;
	int *index2;
public:
	void init(int *_index1, int *_index2) {
		index1 = _index1;
		index2 = _index2;
	}
	bool operator()(int i, int j)
	{
		if (index1[i] < index1[j]) return true;
		if (index1[i] > index1[j]) return false;
		return (index2[i] < index2[j]);
	}
};

void VectorSparseMatrix::sort(SparseMatrix::SortOrder sort_order)
{
printf("VectorSparseMatrix::sort(%d, %ld)\n", sort_order, size());
	// Decide on how we'll sort
	CmpIndex2 cmp;
	switch(sort_order) {
		case SparseMatrix::SortOrder::COLUMN_MAJOR :
			cmp.init(&jndx[0], &indx[0]);
		break;
		default :
			cmp.init(&indx[0], &jndx[0]);
		break;
	}

	// Generate a permuatation
	int n = size();
	std::vector<int> perm; perm.reserve(n);
	for (int i=0; i<n; ++i) perm.push_back(i);
	std::sort(perm.begin(), perm.end(), cmp);

//for (int i=0; i<100; ++i) printf("%d, ", perm[i]); printf("\n");

	// Apply permutation to val
	std::vector<double> dtmp; dtmp.reserve(n);
	for (int i=0; i<n; ++i) dtmp.push_back(val[perm[i]]);
	val = std::move(dtmp);

	// Apply permutation to indx
	std::vector<int> itmp; itmp.reserve(n);
	for (int i=0; i<n; ++i) itmp.push_back(indx[perm[i]]);
	std::swap(itmp, indx);

	// Apply permutation to jndx
	itmp.clear();
	for (int i=0; i<n; ++i) itmp.push_back(jndx[perm[i]]);
	jndx = std::move(itmp);	
}

void VectorSparseMatrix::sum_duplicates(
	SparseMatrix::SortOrder sort_order) // = SparseMatrix::SortOrder::ROW_MAJOR
{
	// Nothing to do for zero-size matrices (or ones with just one element)
	if (size() <= 1) return;

	// Decide on how we'll sort
	CmpIndex2 cmp;
	switch(sort_order) {
		case SparseMatrix::SortOrder::COLUMN_MAJOR :
			cmp.init(&jndx[0], &indx[0]);
		break;
		default :
			cmp.init(&indx[0], &jndx[0]);
		break;
	}

	// Generate a sorted permuatation
	int n = size();
	std::vector<int> perm; perm.reserve(n);
	for (int i=0; i<n; ++i) perm.push_back(i);
	std::sort(perm.begin(), perm.end(), cmp);

	// Output arrays
	std::vector<int> nindx;
	std::vector<int> njndx;
	std::vector<double> nval;

	// Identify duplicates
	nindx.push_back(indx[perm[0]]);
	njndx.push_back(jndx[perm[0]]);
	nval.push_back(val[perm[0]]);
	for (int i=1; i<indx.size(); ++i) {
		int row = indx[perm[i]];
		int col = jndx[perm[i]];
//printf("remove_dup: %d %d\n", row, col);
		if ((row == nindx.back()) && (col == njndx.back())) {
			nval.back() += val[perm[i]];
		} else {
			nindx.push_back(indx[perm[i]]);
			njndx.push_back(jndx[perm[i]]);
			nval.push_back(val[perm[i]]);
		}
	}

	indx = std::move(nindx);
	jndx = std::move(njndx);
	val = std::move(nval);
}


std::unique_ptr<VectorSparseMatrix> VectorSparseMatrix::netcdf_read(
	NcFile &nc, std::string const &vname)
{
	NcVar *indexVar = nc.get_var((vname + ".index").c_str());
	long ne = indexVar->get_dim(0)->size();		// # elements in matrix

	std::vector<int> indx(ne);
	indexVar->set_cur(0,0);
	indexVar->get(&indx[0], ne, 1);

	std::vector<int> jndx(ne);
	indexVar->set_cur(0,1);
	indexVar->get(&jndx[0], ne, 1);

	NcVar *valVar = nc.get_var((vname + ".val").c_str());
	std::vector<double> val(ne);
	valVar->get(&val[0], ne);

	// Read the descriptor
	NcVar *descrVar = nc.get_var((vname + ".descr").c_str());
	SparseDescr descr(
		get_att(descrVar, "nrow")->as_int(0),
		get_att(descrVar, "ncol")->as_int(0),
		get_att(descrVar, "index_base")->as_int(0),
		(SparseMatrix::MatrixStructure)get_att(descrVar, "matrix_structure")->as_int(0),
		(SparseMatrix::TriangularType)get_att(descrVar, "triangular_type")->as_int(0),
		(SparseMatrix::MainDiagonalType)get_att(descrVar, "main_diagonal_type")->as_int(0));

#if 0
// ********* TODO REMOVE ***************
// Fixup bug
for (int i=0; i<ne; ++i) {
	indx[i] -= 1;
	jndx[i] -= 1;
}
#endif

	return std::unique_ptr<VectorSparseMatrix>(new VectorSparseMatrix(descr,
		std::move(indx), std::move(jndx), std::move(val)));
}
// =====================================================================

}	// namespace giss
