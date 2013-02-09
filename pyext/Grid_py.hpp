#pragma once

#include <memory>
#include <glint2/Grid.hpp>

/// Classmembers of the Python class
struct PyGrid {
	PyObject_HEAD
	std::unique_ptr<glint2::Grid> grid;	// The grid, in C++

	// Values copied out of grid
//	long ncells_full;

	void init(std::unique_ptr<glint2::Grid> &&_grid);
//	PyGrid() {}

};

extern PyTypeObject GridType;
