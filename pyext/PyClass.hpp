#pragma once

#include <memory>

/// Classmembers of the Python class
template<class CppT>
struct PyClass {
	PyObject_HEAD
	std::unique_ptr<CppT> ptr;	// The maker, in C++

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
	self->ob_type->tp_free((PyObject *)self);
}

template<class CppT>
PyObject *PyClass<CppT>::new_instance(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	PyClass<CppT> *self;

	self = new (type->tp_alloc(type, 0)) PyClass<CppT>;

//	Py_INCREF(self);
    return (PyObject *)self;
}
