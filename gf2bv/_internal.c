#include "_inernal.h"

// CPython does not have an public API for this yet
// ref: https://github.com/aleaxit/gmpy/issues/467
#if PY_VERSION_HEX >= 0x030C0000
// 3.12 and later
#define NON_SIZE_BITS 3
#define PyLong_DigitCount(o) ((o)->long_value.lv_tag >> NON_SIZE_BITS)
#define GET_OB_DIGITS(o) ((o)->long_value.ob_digit)
#else
// pre 3.12
#define PyLong_DigitCount(o) Py_ABS(Py_SIZE(o))
#define GET_OB_DIGITS(o) ((o)->ob_digit)
#endif

#if PY_VERSION_HEX >= 0x030C0000
// 3.12 and later they are immortal: https://peps.python.org/pep-0683/
#define PythonTrue Py_True
#define PythonFalse Py_False
#else
// pre 3.12
#define PythonTrue Py_NewRef(Py_True)
#define PythonFalse Py_NewRef(Py_False)
#endif

static void nop() {}

static inline PyObject *mzd_vector_to_pylong(char *buf, mzd_t *v) {
	rci_t len = v->ncols;
	// buf[len] must be 0
	for (int c = 0; c < len; c++) {
		buf[len - 1 - c] = mzd_read_bit(v, 0, c) ? '1' : '0';
	}
	return PyLong_FromString(buf, NULL, 2);
}

#define Iter_PyLong_Bits(obj, max, action, after_action)                      \
	do {                                                                      \
		Py_ssize_t _n_digits = PyLong_DigitCount(obj);                        \
		Py_ssize_t bitcnt = 0;                                                \
		for (Py_ssize_t _digits_i = 0; _digits_i < _n_digits && bitcnt < max; \
		     _digits_i++) {                                                   \
			digit _digit = GET_OB_DIGITS(obj)[_digits_i];                     \
			for (int _digit_j = 0; _digit_j < PyLong_SHIFT && bitcnt < max;   \
			     _digit_j++) {                                                \
				int bit = _digit & 1;                                         \
				action;                                                       \
				_digit >>= 1;                                                 \
				bitcnt++;                                                     \
			}                                                                 \
		}                                                                     \
		for (; bitcnt < max; bitcnt++) {                                      \
			after_action;                                                     \
		}                                                                     \
	} while (0)

#pragma region AffineSpaceIterable

PyObject *affinespaceiterable_next(AffineSpaceIterableObject *self) {
	// TODO: use gray code for more efficient enumeration

	mzd_t *result;
	rci_t n, r;
	int sentinel;
	AffineSpaceObject *space = self->space;

	n = space->basis->nrows;
	if (self->state[n])
		return NULL; /* StopIteration */

	result = mzd_copy(NULL, space->origin);
	for (r = 0; r < n; r++)
		if (self->state[r])
			mzd_combine_even_in_place(result, 0, 0, space->basis, r, 0);

	sentinel = 1;
	for (r = 0; r < n; r++) {
		if ((self->state[r] ^= 1)) {
			sentinel = 0;
			break;
		}
	}
	self->state[n] = sentinel;

	PyObject *ret = mzd_vector_to_pylong(self->str, result);
	mzd_free(result);
	return ret;
}

static void affinespaceiterable_dealloc(AffineSpaceIterableObject *self) {
	PyObject_GC_UnTrack(self);
	Py_DECREF(self->space);
	free(self->state);
	free(self->str);
	PyObject_GC_Del(self);
}

static int affinespaceiterable_traverse(AffineSpaceIterableObject *self,
                                        visitproc visit,
                                        void *arg) {
	Py_VISIT(self->space);
	return 0;
}

static int affinespaceiterable_clear(AffineSpaceIterableObject *self) {
	Py_CLEAR(self->space);
	return 0;
}

static PyTypeObject AffineSpaceIterable_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name =
        "_internal.AffineSpaceIterable",
    .tp_basicsize = sizeof(AffineSpaceIterableObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)affinespaceiterable_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)affinespaceiterable_traverse,
    .tp_clear = (inquiry)affinespaceiterable_clear,
    .tp_doc = NULL,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)affinespaceiterable_next,
};

#pragma endregion

#pragma region AffineSpace

static PyObject *affinespace_iter(PyObject *self) {
	// create an iterator object
	AffineSpaceIterableObject *it =
	    PyObject_GC_New(AffineSpaceIterableObject, &AffineSpaceIterable_Type);
	it->space = (AffineSpaceObject *)self;
	it->state = malloc(it->space->basis->nrows + 1);
	it->str = malloc(it->space->origin->ncols + 1);
	it->str[it->space->origin->ncols] = '\0';
	memset(it->state, 0, it->space->basis->nrows + 1);
	Py_INCREF(self);
	PyObject_GC_Track(it);
	return (PyObject *)it;
}

static PyObject *affinespace_get_dimension(AffineSpaceObject *self) {
	return PyLong_FromLong(self->basis->nrows);
}

static PyObject *affinespace_get_origin(AffineSpaceObject *self) {
	char *str = malloc(self->origin->ncols + 1);
	str[self->origin->ncols] = '\0';
	PyObject *ret = mzd_vector_to_pylong(str, self->origin);
	free(str);
	return ret;
}

static PyObject *affinespace_get_basis(AffineSpaceObject *self) {
	PyObject *ret = PyTuple_New(self->basis->nrows);
	char *str = malloc(self->basis->ncols + 1);
	for (rci_t r = 0; r < self->basis->nrows; r++) {
		str[self->basis->ncols] = '\0';
		mzd_t *w =
		    mzd_init_window(self->basis, r, 0, r + 1, self->basis->ncols);
		PyTuple_SET_ITEM(ret, r, mzd_vector_to_pylong(str, w));
		mzd_free_window(w);
	}
	free(str);
	return ret;
}

static PyGetSetDef AffineSpace_getsetters[] = {
    {"dimension", (getter)affinespace_get_dimension, NULL,
     "Dimension of the affine space", NULL},
    {"origin", (getter)affinespace_get_origin, NULL,
     "Origin of the affine space", NULL},
    {"basis", (getter)affinespace_get_basis, NULL, "Basis of the affine space",
     NULL},
    {NULL} /* Sentinel */
};

static PyObject *affinespace_get(AffineSpaceObject *self,
                                 PyObject *const *args,
                                 Py_ssize_t nargs) {
	PyObject *index;
	if (nargs != 1) {
		PyErr_SetString(PyExc_TypeError, "get requires 1 argument");
	}
	index = args[0];
	if (!PyLong_Check(index)) {
		PyErr_SetString(PyExc_TypeError, "Index must be an integer");
		return NULL;
	}
	mzd_t *result = mzd_copy(NULL, self->origin);
	rci_t nr = self->basis->nrows;

	PyLongObject *n = (PyLongObject *)index;
	Iter_PyLong_Bits(n, nr,
	                 {
		                 if (bit) {
			                 mzd_combine_even_in_place(result, 0, 0,
			                                           self->basis, bitcnt, 0);
		                 }
	                 },
	                 {});

	char *str = malloc(result->ncols + 1);
	str[result->ncols] = '\0';
	PyObject *ret = mzd_vector_to_pylong(str, result);
	free(str);
	mzd_free(result);
	return ret;
}

static PyMethodDef AffineSpace_methods[] = {
    {"get", _PyCFunction_CAST(affinespace_get), METH_FASTCALL,
     "get(n)\n"
     "--\n"
     "\n"
     "Get the n-th element of the affine space, should check 0 <= n < "
     "2**(space.dimension) first."},
    {NULL} /* Sentinel */
};

static void affinespace_dealloc(AffineSpaceObject *self) {
	PyObject_GC_UnTrack(self);
	mzd_free(self->origin);
	mzd_free(self->basis);
	PyObject_GC_Del(self);
}

static PyTypeObject AffineSpace_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name = "_internal.AffineSpace",
    .tp_basicsize = sizeof(AffineSpaceObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)affinespace_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)nop,  // no Python objects to traverse
    .tp_clear = (inquiry)nop,          // no Python objects to clear
    .tp_doc = NULL,
    .tp_iter = affinespace_iter,
    .tp_methods = AffineSpace_methods,
    .tp_getset = AffineSpace_getsetters,
};

#pragma endregion

PyObject *m4ri_solve(PyObject *self, PyObject *const *args, Py_ssize_t nargs) {
	PyObject *linsys_list;
	Py_ssize_t cols;
	long mode = 0;
	if (nargs != 3) {
		PyErr_SetString(PyExc_TypeError, "m4ri_solve requires 3 arguments");
		return NULL;
	}
	linsys_list = args[0];
	if (!PyList_Check(linsys_list)) {
		PyErr_SetString(PyExc_TypeError,
		                "The first argument equations must be a list");
		return NULL;
	}
	cols = PyLong_AsSsize_t(args[1]);
	if (cols <= 0) {
		if (cols == -1 && PyErr_Occurred()) {
			return NULL;
		}
		PyErr_SetString(PyExc_ValueError, "Number of columns must be positive");
		return NULL;
	}
	mode = PyLong_AsLong(args[2]);
	if (mode == -1 && PyErr_Occurred()) {
		return NULL;
	}
	if (mode != SOLVE_MODE_SINGLE && mode != SOLVE_MODE_ALL) {
		PyErr_SetString(PyExc_ValueError, "Invalid mode");
		return NULL;
	}
	int need_kernel = mode == SOLVE_MODE_ALL;
	Py_ssize_t rows = PyList_Size(linsys_list);
	if (rows < cols) {
		PyErr_SetString(PyExc_ValueError,
		                "Number of rows must be greater than or equal to "
		                "number of columns, try pad with zeros.");
		return NULL;
	}

	// Ax = B
	mzd_t *A = mzd_init(rows, cols);
	mzd_t *B = mzd_init(rows, 1);

	int B_is_not_zero = 0;

	for (Py_ssize_t r = 0; r < rows; r++) {
		PyObject *item = PyList_GetItem(linsys_list, r);
		if (!PyLong_Check(item)) {
			mzd_free(A);
			mzd_free(B);
			PyErr_SetString(PyExc_TypeError, "List items must be integers");
			return NULL;
		}
		// each item's lsb is the affine term
		// and higher bits is the linear terms (total: cols)
		// row in A: [v_1 v_2 ... v_cols]
		// row in B: [v_0]
		PyLongObject *v = (PyLongObject *)item;
		Iter_PyLong_Bits(v, cols + 1,
		                 {
			                 if (bitcnt == 0) {
				                 mzd_write_bit(B, r, 0, bit);
				                 B_is_not_zero |= bit;
			                 } else {
				                 mzd_write_bit(A, r, bitcnt - 1, bit);
			                 }
		                 },
		                 {});
	}

	// first, find the base solution
	mzd_t *origin;
	{
		mzd_t *A_copy = 0;
		if (B_is_not_zero) {
			if (need_kernel) {
				A_copy = mzd_copy(NULL, A);
			}
			// printf("A->nrows: %d, A->ncols: %d\n", A->nrows, A->ncols);
			// printf("B->nrows: %d, B->ncols: %d\n", B->nrows, B->ncols);
			// this is actually solve_right in Sage...
			if (mzd_solve_left(A, B, 0, 1) != 0) {
				mzd_free(A);
				mzd_free(B);
				if (need_kernel) {
					mzd_free(A_copy);
				}
				return Py_None;
			}
			// the base solution is stored in B as column vector
			B->nrows = cols;
			origin = mzd_transpose(0, B);
			mzd_free(A);  // A is overwritten by mzd_solve_left
			A = A_copy;
		} else {
			// the trivial solution
			origin = mzd_init(1, cols);
		}
	}
	mzd_free(B);
	// if need_kernel: A, origin is valid
	// else: A is valid or 0, origin is valid

	// if we only need one solution, return it
	if (mode == SOLVE_MODE_SINGLE) {
		char *str = malloc(cols + 1);
		str[cols] = '\0';
		PyObject *ret = mzd_vector_to_pylong(str, origin);
		free(str);
		mzd_free(origin);
		if (A) {
			mzd_free(A);
		}
		return ret;
	}

	// compute the basis
	mzd_t *tker;
	{
		// and this is right_kernel_matrix
		mzd_t *ker = mzd_kernel_left_pluq(A, 0);
		mzd_free(A);
		if (ker == NULL) {
			tker = mzd_init(0, 0);
		} else {
			tker = mzd_transpose(0, ker);
			mzd_free(ker);
		}
	}
	// now, only origin and tker is valid

	AffineSpaceObject *it =
	    PyObject_GC_New(AffineSpaceObject, &AffineSpace_Type);
	it->origin = origin;
	it->basis = tker;
	PyObject_GC_Track(it);
	return (PyObject *)it;
}

PyObject *to_bits(PyObject *self, PyObject *const *args, Py_ssize_t nargs) {
	Py_ssize_t n;
	PyObject *a;
	if (nargs != 2) {
		PyErr_SetString(PyExc_TypeError, "to_bits requires 2 arguments");
		return NULL;
	}
	n = PyLong_AsSsize_t(args[0]);
	if (n < 0) {
		if (n == -1 && PyErr_Occurred()) {
			return NULL;
		}
		PyErr_SetString(PyExc_ValueError, "n must be non-negative");
		return NULL;
	}
	a = args[1];
	if (!PyLong_Check(a)) {
		PyErr_SetString(PyExc_TypeError, "a must be an integer");
		return NULL;
	}
	PyObject *list = PyList_New(n);
	PyLongObject *a_long = (PyLongObject *)a;
	Iter_PyLong_Bits(
	    a_long, n,
	    { PyList_SetItem(list, bitcnt, bit ? PythonTrue : PythonFalse); },
	    { PyList_SetItem(list, bitcnt, PythonFalse); });
	return list;
}

static void set_bits(char *bits, Py_ssize_t n, PyLongObject *o) {
	Iter_PyLong_Bits(
	    o, n, { bits[bitcnt] = bit ? 1 : 0; }, { bits[bitcnt] = 0; });
}

PyObject *mul_bit_quad(PyObject *self,
                       PyObject *const *args,
                       Py_ssize_t nargs) {
	Py_ssize_t n;
	PyObject *a, *b, *v, *basis;
	// if (!PyArg_ParseTuple(args, "nO!O!O!O!", &n, &PyLong_Type, &a, &PyLong_Type,
	//                       &b, &PyLong_Type, &v, &PyList_Type, &basis))
	// 	return NULL;
	if (nargs != 5) {
		PyErr_SetString(PyExc_TypeError, "mul_bit_quad requires 5 arguments");
		return NULL;
	}
	n = PyLong_AsSsize_t(args[0]);
	if (n <= 0) {
		if (n == -1 && PyErr_Occurred()) {
			return NULL;
		}
		PyErr_SetString(PyExc_ValueError, "n must be positive");
		return NULL;
	}
	a = args[1];
	b = args[2];
	v = args[3];
	if (!PyLong_Check(a) || !PyLong_Check(b) || !PyLong_Check(v)) {
		PyErr_SetString(PyExc_TypeError, "a and b and v must be integers");
		return NULL;
	}
	basis = args[4];
	if (!PyList_Check(basis)) {
		PyErr_SetString(PyExc_TypeError, "basis must be a list");
		return NULL;
	}
	Py_ssize_t len_basis = PyList_GET_SIZE(basis);
	if (len_basis != 1 + n + n * (n - 1) / 2) {
		PyErr_SetString(PyExc_ValueError, "The length of basis is not correct");
		return NULL;
	}
	char *a_bits = malloc(n);
	set_bits(a_bits, n, (PyLongObject *)a);
	char *b_bits = malloc(n);
	set_bits(b_bits, n, (PyLongObject *)b);

	Py_INCREF(v);  // invariant: v is a strong reference

	int mi = 1 + n;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			int r = (a_bits[i] & b_bits[j]) ^ (a_bits[j] & b_bits[i]);
			if (r) {
				PyObject *ret = PyNumber_Or(v, PyList_GetItem(basis, mi));
				if (ret == NULL) {
					free(a_bits);
					free(b_bits);
					PyErr_SetString(PyExc_TypeError,
					                "Failed to compute or, list items "
					                "must be integers");
					return NULL;
				}
				Py_SETREF(v, ret);
			}
			mi++;
		}
	}
	free(a_bits);
	free(b_bits);
	return v;
}

static PyMethodDef methods[] = {
    {"m4ri_solve", _PyCFunction_CAST(m4ri_solve), METH_FASTCALL,
     "m4ri_solve(equations, cols, mode)\n"
     "--\n"
     "\n"
     "Solve a linear system over GF(2) with M4RI"},
    {"to_bits", _PyCFunction_CAST(to_bits), METH_FASTCALL,
     "to_bits(n, number)\n"
     "--\n"
     "\n"
     "Convert an integer to a list of bits (bool values)"},
    {"mul_bit_quad", _PyCFunction_CAST(mul_bit_quad), METH_FASTCALL,
     "mul_bit_quad(n, a, b, v, basis)\n"
     "--\n"
     "\n"
     "Multiply two linear symbolic bits to a linearized quadratic symbolic "
     "bit"},
    {NULL} /* Sentinel */
};

static struct PyModuleDef _internal = {PyModuleDef_HEAD_INIT, "_internal", NULL,
                                       -1, methods};

#define INIT_TYPE(type)              \
	do {                             \
		if (PyType_Ready(&type) < 0) \
			return NULL;             \
	} while (0)

#define ADD_TYPE(mod, type)                     \
	do {                                        \
		if (PyModule_AddType(mod, &type) < 0) { \
			Py_DECREF(mod);                     \
			return NULL;                        \
		}                                       \
	} while (0)

PyMODINIT_FUNC PyInit__internal(void) {
	INIT_TYPE(AffineSpace_Type);
	INIT_TYPE(AffineSpaceIterable_Type);
	PyObject *mod = PyModule_Create(&_internal);
	if (mod == NULL)
		return NULL;
	ADD_TYPE(mod, AffineSpace_Type);
	ADD_TYPE(mod, AffineSpaceIterable_Type);
	return mod;
}
