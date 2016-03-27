/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <memory>

/// Classmembers of the Python class
template<class CppT>
struct PyClass {
    PyObject_HEAD
    std::unique_ptr<CppT> ptr;  // The maker, in C++

    void init(std::unique_ptr<CppT> &&_ptr);

    // For use in Python type descriptor
    static PyObject *new_instance(PyTypeObject *type, PyObject *args, PyObject *kwds);
    static void dealloc(PyClass<CppT> *self);

};

template<class CppT>
void PyClass<CppT>::init(std::unique_ptr<CppT> &&_ptr)
{
    ptr = std::move(_ptr);
}

//template<class CppT>
//void PyClass_dealloc(PyClass<CppT> *self);

template<class CppT>
void PyClass<CppT>::dealloc(PyClass<CppT> *self)
{
    self->~PyClass<CppT>();
    self->ob_base.ob_type->tp_free((PyObject *)self);
}

template<class CppT>
PyObject *PyClass<CppT>::new_instance(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyClass<CppT> *self;

    self = new (type->tp_alloc(type, 0)) PyClass<CppT>;

//  Py_INCREF(self);
    return (PyObject *)self;
}
