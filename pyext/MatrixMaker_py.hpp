#pragma once

#include <memory>
#include <glint2/MatrixMaker.hpp>

/// Classmembers of the Python class
struct PyMatrixMaker {
	PyObject_HEAD
	std::unique_ptr<glint2::MatrixMaker> maker;	// The maker, in C++

	void init(std::unique_ptr<glint2::MatrixMaker> &&_maker);

};

extern PyTypeObject MatrixMakerType;
