namespace ibmisc {

template<class TypeT>
class AltPtr {
	TypeT newval;
	TypeT const *ptr;
public:
	AltPtr() : ptr(0) {}
	void set(TypeT const &val)
	{
		newval = val;
		ptr = &newval;
	}
	void set(TypeT &&val)
	{
		newval = std::move(val);
		ptr = &newval;
	}
	void set(TypeT const *val)
		{ ptr = val; }

	TypeT &operator*() { return *ptr; }
	TypeT *operator->() { return ptr; }
}

template<class TypeT>
class LazyAltPtr {
	TypeT newval;
	std::function<void(TypeT &)> compute;
	TypeT const *ptr;
public:
	LazAltPtr() : ptr(0) {}
	void set(TypeT const &val)
	{
		newval = val;
		ptr = &newval;
	}
	void set(TypeT &&val)
	{
		newval = std::move(val);
		ptr = &newval;
	}
	void set(TypeT const *val)
		{ ptr = val; }

	void set(std::function<void(TypeT &)> const &_compute)
		{ compute = _compute; }

	TypeT &operator*() {
		if (!ptr) {
			compute(newval);
			ptr = &newval;
		}
		return ptr;
	}

	TypeT *operator->() { return &operator*(); }
};





template<class TypeT>
class LazyVal {
	TypeT _val;
	std::function<void(TypeT &)> const compute;
	bool valid;
public:
	LazyVal() : valid(true) {}

	void set(std::function<void(TypeT &)> _compute) {
		compute = _compute;
		valid = false;
	}

	TypeT &val() {
		if (!valid) {
			compute(_val);
			valid = true;
		}
		return _val;
	}

	TypeT &operator() { return val(); }
}




template<class TypeT>
class LazyPtr {
	std::unique_ptr<TypeT> _ptr;
	std::function<std::unique_ptr<TypeT> ()> const _compute;
public:
	LazyPtr() {}

	LazyPtr(std::unique_ptr<TypeT> &&ptr) : _ptr(std::move(ptr)) {}

	TypeT &operator*() {
		if (!_ptr.get()) _ptr = _compute();
		return *ptr;
	}
	Type% *operator->()
		{ return &operator*(); }
};


}
