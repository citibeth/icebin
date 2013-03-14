#pragma once

namespace giss {

template<class T>
class DynArray {
public:
	size_t const size;
	size_t const ele_size;
	char * const buf;
	char * const buf_end;

	DynArray(size_t _ele_size, size_t _size) :
		ele_size(_ele_size),
		size(_size),
		buf((char *)malloc(size * ele_size)),
		buf_end(buf + size * ele_size) {}

	~DynArray() { free(buf); }

	DynArray(DynArray const &b) = delete;
	void operator=(DynArray const &b) = delete;
	void operator=(DynArray &&b) = delete;

	T *begin() { return reinterpret_cast<T *>(buf); }
	T *end() { return  reinterpret_cast<T *>(buf_end); }

	T &operator[](int ix)
		{ return *reinterpret_cast<T *>(buf + ele_size * ix); }

	T const &operator[](int ix) const
		{ return *reinterpret_cast<T const *>(buf + ele_size * ix); }

	T *incr(T *ptr) const {
		return reinterpret_cast<T *>(
			reinterpret_cast<char *>(ptr) + ele_size);
	}

	T *decr(T *ptr) const {
		return reinterpret_cast<T *>(
			reinterpret_cast<char *>(ptr) - ele_size);
	}

};

}
