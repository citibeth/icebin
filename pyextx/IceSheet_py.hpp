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

#pragma once

#include <memory>
#include <glint2/IceSheet.hpp>

/// Classmembers of the Python class
struct PyIceSheet {
    PyObject_HEAD
    std::unique_ptr<glint2::IceSheet> sheet;    // The sheet, in C++

    void init(std::unique_ptr<glint2::IceSheet> &&_sheet);

};

extern PyTypeObject IceSheetType;
