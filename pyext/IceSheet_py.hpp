#pragma once

#include <memory>
#include <glint2/IceSheet.hpp>

/// Classmembers of the Python class
struct PyIceSheet {
	PyObject_HEAD
	std::unique_ptr<glint2::IceSheet> sheet;	// The sheet, in C++

	void init(std::unique_ptr<glint2::IceSheet> &&_sheet);

};

extern PyTypeObject IceSheetType;
