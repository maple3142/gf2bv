#include "_internal.h"

#include <dlfcn.h>

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
#define PythonNone Py_None
#else
// pre 3.12
#define PythonTrue Py_NewRef(Py_True)
#define PythonFalse Py_NewRef(Py_False)
#define PythonNone Py_NewRef(Py_None)
#endif

static PyTypeObject GF2Matrix_Type;

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

#pragma region AffineSpaceIterator

static PyObject *affinespaceiterator_next_naive(
    AffineSpaceIteratorObject *self) {
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

static void affinespaceiterator_dealloc_naive(AffineSpaceIteratorObject *self) {
	PyObject_GC_UnTrack(self);
	Py_DECREF(self->space);
	free(self->state);
	free(self->str);
	PyObject_GC_Del(self);
}

static PyObject *affinespaceiterator_next_graycode(
    AffineSpaceIteratorObject *self) {
	if (self->gray.cur == NULL) {
		return NULL; /* StopIteration */
	}
	PyObject *ret = mzd_vector_to_pylong(self->str, self->gray.cur);

	// use gray code to compute the next vector
	uint64_t x = (self->gray.idx) ^ (self->gray.idx >> 1);
	self->gray.idx++;
	uint64_t y = (self->gray.idx) ^ (self->gray.idx >> 1);
	int diff_idx = __builtin_ctzll(x ^ y);
	if (diff_idx >= self->space->basis->nrows ||
	    (self->space->basis->nrows == 64 && self->gray.idx == 0)) {
		mzd_free(self->gray.cur);
		self->gray.cur = NULL;
		return ret;
	}
	mzd_combine_even_in_place(self->gray.cur, 0, 0, self->space->basis,
	                          diff_idx, 0);
	return ret;
}

static void affinespaceiterator_dealloc_graycode(
    AffineSpaceIteratorObject *self) {
	PyObject_GC_UnTrack(self);
	Py_DECREF(self->space);
	if (self->gray.cur)
		// if enumeration ended, self->gray.cur is NULL
		mzd_free(self->gray.cur);
	free(self->str);
	PyObject_GC_Del(self);
}

static int affinespaceiterator_traverse(AffineSpaceIteratorObject *self,
                                        visitproc visit,
                                        void *arg) {
	Py_VISIT(self->space);
	return 0;
}

static int affinespaceiterator_clear(AffineSpaceIteratorObject *self) {
	Py_CLEAR(self->space);
	return 0;
}

// only differes in tp_iternext and tp_dealloc

static PyTypeObject AffineSpaceIterator_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name =
        "_internal.AffineSpaceIterator",
    .tp_basicsize = sizeof(AffineSpaceIteratorObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)affinespaceiterator_dealloc_graycode,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)affinespaceiterator_traverse,
    .tp_clear = (inquiry)affinespaceiterator_clear,
    .tp_doc = NULL,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)affinespaceiterator_next_graycode,
};

static PyTypeObject AffineSpaceIteratorSlow_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name =
        "_internal.AffineSpaceIteratorSlow",
    .tp_basicsize = sizeof(AffineSpaceIteratorObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)affinespaceiterator_dealloc_naive,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)affinespaceiterator_traverse,
    .tp_clear = (inquiry)affinespaceiterator_clear,
    .tp_doc = NULL,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)affinespaceiterator_next_naive,
};

#pragma endregion

#pragma region AffineSpace

static PyObject *affinespace_iter(PyObject *self) {
	AffineSpaceObject *space = (AffineSpaceObject *)self;
	char *str = malloc(space->origin->ncols + 1);
	str[space->origin->ncols] = '\0';
	int use_gray = space->basis->nrows <= 64;
	PyTypeObject *type =
	    use_gray ? &AffineSpaceIterator_Type : &AffineSpaceIteratorSlow_Type;

	// create an iterator object
	AffineSpaceIteratorObject *it =
	    PyObject_GC_New(AffineSpaceIteratorObject, type);
	it->space = space;
	it->str = str;
	if (use_gray) {
		it->gray.cur = mzd_copy(NULL, it->space->origin);
		it->gray.idx = 0;
	} else {
		it->state = malloc(space->basis->nrows + 1);
		memset(it->state, 0, space->basis->nrows + 1);
	}
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

// https://github.com/malb/m4ri/blob/775189bfea96ffaeab460513413fcf4fbcd64392/m4ri/solve.c#L148
mzd_t *_mzd_kernel_left_pluq(mzd_t *A,
                             rci_t r,
                             mzp_t *P,
                             mzp_t *Q,
                             int const cutoff) {
	//   mzp_t *P = mzp_init(A->nrows);
	//   mzp_t *Q = mzp_init(A->ncols);

	//   rci_t r = mzd_pluq(A, P, Q, cutoff);

	//   if (r == A->ncols) {
	//     mzp_free(P);
	//     mzp_free(Q);
	//     __M4RI_DD_MZD(A);
	//     return NULL;
	//   }

	if (r == A->ncols) {
		return NULL;
	}

	mzd_t *U = mzd_init_window(A, 0, 0, r, r);

	mzd_t *R = mzd_init(A->ncols, A->ncols - r);
	mzd_t *RU = mzd_init_window(R, 0, 0, r, R->ncols);

	for (rci_t i = 0; i < r; i++) {
		for (rci_t j = 0; j < RU->ncols; j += m4ri_radix) {
			const int workload = MIN(m4ri_radix, RU->ncols - j);
			mzd_xor_bits(RU, i, j, workload,
			             mzd_read_bits(A, i, r + j, workload));
		}
	}

	mzd_trsm_upper_left(U, RU, cutoff);

	for (rci_t i = 0; i < R->ncols; ++i) {
		mzd_write_bit(R, r + i, i, 1);
	}
	mzd_apply_p_left_trans(R, Q);
	// mzp_free(P);
	// mzp_free(Q);
	mzd_free_window(RU);
	mzd_free_window(U);

	// __M4RI_DD_MZD(A);
	// __M4RI_DD_MZD(R);
	return R;
}

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
	if (mode != SOLVE_MODE_SINGLE && mode != SOLVE_MODE_AFFINE_SPACE) {
		PyErr_SetString(PyExc_ValueError, "Invalid mode");
		return NULL;
	}
	Py_ssize_t rows = PyList_GET_SIZE(linsys_list);
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
		PyObject *item = PyList_GET_ITEM(linsys_list, r);
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

	// release GIL because we don't need to call Python API now
	PyThreadState *_save = PyEval_SaveThread();

	mzp_t *P = mzp_init(A->nrows);
	mzp_t *Q = mzp_init(A->ncols);
	rci_t r = _mzd_pluq(A, P, Q, 0);

	// first, find the base solution
	mzd_t *sol0;
	{
		if (B_is_not_zero) {
			// this is actually solve_right in Sage...
			if (_mzd_pluq_solve_left(A, r, P, Q, B, 0, 1) != 0) {
				mzd_free(A);
				mzd_free(B);
				mzp_free(P);
				mzp_free(Q);
				PyEval_RestoreThread(_save);
				return PythonNone;
			}
			// the base solution is stored in B as column vector
			B->nrows = cols;
			sol0 = mzd_transpose(0, B);
		} else {
			// the trivial solution
			sol0 = mzd_init(1, cols);
		}
	}
	mzd_free(B);
	// A, sol0, P, Q is valid

	// if we only need one solution, return it
	if (mode == SOLVE_MODE_SINGLE) {
		PyEval_RestoreThread(_save);
		char *str = malloc(cols + 1);
		str[cols] = '\0';
		PyObject *ret = mzd_vector_to_pylong(str, sol0);
		free(str);

		mzd_free(A);
		mzd_free(sol0);
		mzp_free(P);
		mzp_free(Q);
		return ret;
	}

	// compute the basis
	mzd_t *tker;
	{
		// and this is right_kernel_matrix
		// mzd_t *ker = mzd_kernel_left_pluq(A, 0);
		mzd_t *ker = _mzd_kernel_left_pluq(A, r, P, Q, 0);
		mzd_free(A);
		mzp_free(P);
		mzp_free(Q);
		if (ker == NULL) {
			tker = mzd_init(0, 0);
		} else {
			tker = mzd_transpose(0, ker);
			mzd_free(ker);
		}
	}
	// now, only sol0 and tker is valid

	PyEval_RestoreThread(_save);

	// create the affine space object and return it
	// sol0 and tker are tracked by it
	AffineSpaceObject *it =
	    PyObject_GC_New(AffineSpaceObject, &AffineSpace_Type);
	it->origin = sol0;
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

PyObject *xor_list(PyObject *self, PyObject *const *args, Py_ssize_t nargs) {
	PyObject *a, *b;
	if (nargs != 2) {
		PyErr_SetString(PyExc_TypeError, "xor_list requires 2 arguments");
		return NULL;
	}
	a = args[0];
	b = args[1];
	if (!PyList_Check(a) || !PyList_Check(b)) {
		PyErr_SetString(PyExc_TypeError, "a and b must be lists");
		return NULL;
	}
	Py_ssize_t len_a = PyList_GET_SIZE(a);
	Py_ssize_t len_b = PyList_GET_SIZE(b);
	if (len_a != len_b) {
		PyErr_SetString(PyExc_ValueError, "The length of a and b is not equal");
		return NULL;
	}
	PyObject *ret = PyList_New(len_a);
	for (Py_ssize_t i = 0; i < len_a; i++) {
		PyObject *item_a = PyList_GET_ITEM(a, i);
		PyObject *item_b = PyList_GET_ITEM(b, i);
		PyObject *xor_item = PyNumber_Xor(item_a, item_b);
		if (xor_item == NULL) {
			Py_DECREF(ret);
			PyErr_SetString(
			    PyExc_TypeError,
			    "Failed to compute xor, list items must be integers");
			return NULL;
		}
		PyList_SET_ITEM(ret, i, xor_item);
	}
	return ret;
}

PyObject *list_where(PyObject *self, PyObject *const *args, Py_ssize_t nargs) {
	PyObject *cond, *a, *b;
	if (nargs != 3) {
		PyErr_SetString(PyExc_TypeError, "list_where requires 3 arguments");
		return NULL;
	}
	cond = args[0];
	if (!PyList_Check(cond)) {
		PyErr_SetString(PyExc_TypeError, "cond must be a list");
		return NULL;
	}
	a = args[1];
	b = args[2];
	int a_is_list = PyList_Check(a);
	int b_is_list = PyList_Check(b);
	Py_ssize_t len_cond = PyList_GET_SIZE(cond);
	if (a_is_list && PyList_GET_SIZE(a) != len_cond) {
		PyErr_SetString(PyExc_ValueError,
		                "The length of a and cond is not equal");
		return NULL;
	}
	if (b_is_list && PyList_GET_SIZE(b) != len_cond) {
		PyErr_SetString(PyExc_ValueError,
		                "The length of b and cond is not equal");
		return NULL;
	}
	for (Py_ssize_t i = 0; i < len_cond; i++) {
		PyObject *item_cond = PyList_GET_ITEM(cond, i);
		PyObject *item_a = a_is_list ? PyList_GET_ITEM(a, i) : a;
		PyObject *item_b = b_is_list ? PyList_GET_ITEM(b, i) : b;
		PyObject *ret_item = PyObject_IsTrue(item_cond) ? item_a : item_b;
		PyList_SET_ITEM(cond, i, Py_NewRef(ret_item));
	}
	Py_IncRef(cond);
	return cond;
}

// gd.h functions needed by eqs_to_sage_mat_helper
void *(*gdImageCreate)(int, int);
int (*gdImageColorAllocate)(void *, int, int, int);
void (*gdImageFilledRectangle)(void *, int, int, int, int, int);
void (*gdImageSetPixel)(void *, int, int, int);
void *(*gdImagePngPtrEx)(void *, int *, int);
void (*gdFree)(void *);
void (*gdImageDestroy)(void *);

PyObject *eqs_to_sage_mat_helper(PyObject *self,
                                 PyObject *const *args,
                                 Py_ssize_t nargs) {
	if (!gdImageCreate) {
		// load gd library dynamically
		void *gd = dlopen("libgd.so", RTLD_LAZY);
		if (!gd) {
			PyErr_Format(PyExc_RuntimeError, "Failed to load libgd.so: %s",
			             dlerror());
			return NULL;
		}
		gdImageCreate = dlsym(gd, "gdImageCreate");
		gdImageColorAllocate = dlsym(gd, "gdImageColorAllocate");
		gdImageFilledRectangle = dlsym(gd, "gdImageFilledRectangle");
		gdImageSetPixel = dlsym(gd, "gdImageSetPixel");
		gdImagePngPtrEx = dlsym(gd, "gdImagePngPtrEx");
		gdFree = dlsym(gd, "gdFree");
		gdImageDestroy = dlsym(gd, "gdImageDestroy");
		if (!gdImageCreate || !gdImageColorAllocate ||
		    !gdImageFilledRectangle || !gdImageSetPixel || !gdImagePngPtrEx ||
		    !gdFree || !gdImageDestroy) {
			PyErr_SetString(PyExc_RuntimeError,
			                "Failed to load functions from libgd.so");
			return NULL;
		}
	}
	// see https://github.com/sagemath/sage/blob/7726cd9e1d01ad32b0f14374c9a4096989c87e14/src/sage/matrix/matrix_mod2_dense.pyx#L1764-L1808
	PyObject *linsys_list;
	Py_ssize_t cols;
	if (nargs != 2) {
		PyErr_SetString(PyExc_TypeError,
		                "eqs_to_sage_mat_helper requires 2 arguments");
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
	Py_ssize_t rows = PyList_GET_SIZE(linsys_list);
	PyObject *affine = PyList_New(rows);

	void *im = gdImageCreate(cols, rows);
	int black = gdImageColorAllocate(im, 0, 0, 0);
	int white = gdImageColorAllocate(im, 255, 255, 255);
	gdImageFilledRectangle(im, 0, 0, cols - 1, rows - 1, white);
	for (Py_ssize_t i = 0; i < rows; i++) {
		PyLongObject *v = (PyLongObject *)PyList_GET_ITEM(linsys_list, i);
		Iter_PyLong_Bits(v, cols + 1,
		                 {
			                 if (bitcnt == 0) {
				                 PyList_SET_ITEM(
				                     affine, i, bit ? PythonTrue : PythonFalse);
			                 }
			                 if (bit) {
				                 gdImageSetPixel(im, bitcnt - 1, i, black);
			                 }
		                 },
		                 {});
	}
	int size;
	char *buf = gdImagePngPtrEx(im, &size, 0);
	PyObject *png_bytes = PyBytes_FromStringAndSize(buf, size);
	gdFree(buf);
	gdImageDestroy(im);
	PyObject *ret = PyTuple_Pack(2, png_bytes, affine);
	Py_DECREF(affine);
	Py_DECREF(png_bytes);
	return ret;
}

PyObject *eqs_to_linear_system(PyObject *self,
                               PyObject *const *args,
                               Py_ssize_t nargs) {
	PyObject *linsys_list;
	Py_ssize_t cols;
	if (nargs != 2) {
		PyErr_SetString(PyExc_TypeError,
		                "eqs_to_linear_system requires 2 arguments");
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
	Py_ssize_t rows = PyList_GET_SIZE(linsys_list);

	mzd_t *A = mzd_init(rows, cols);
	mzd_t *B = mzd_init(rows, 1);

	for (Py_ssize_t r = 0; r < rows; r++) {
		PyObject *item = PyList_GET_ITEM(linsys_list, r);
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
			                 } else {
				                 mzd_write_bit(A, r, bitcnt - 1, bit);
			                 }
		                 },
		                 {});
	}
	GF2MatrixObject *Aobj = PyObject_GC_New(GF2MatrixObject, &GF2Matrix_Type);
	Aobj->matrix = A;
	GF2MatrixObject *Bobj = PyObject_GC_New(GF2MatrixObject, &GF2Matrix_Type);
	Bobj->matrix = B;
	PyObject_GC_Track(Aobj);
	PyObject_GC_Track(Bobj);
	return PyTuple_Pack(2, (PyObject *)Aobj, (PyObject *)Bobj);
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
    {"xor_list", _PyCFunction_CAST(xor_list), METH_FASTCALL,
     "xor_list(a, b)\n"
     "--\n"
     "\n"
     "XOR two lists of integers"},
    {"list_where", _PyCFunction_CAST(list_where), METH_FASTCALL,
     "list_where(cond, a, b)\n"
     "--\n"
     "\n"
     "Select elements from a or b based on the condition like np.where, and "
     "set it to the cond "
     "list"},
    {"eqs_to_sage_mat_helper", _PyCFunction_CAST(eqs_to_sage_mat_helper),
     METH_FASTCALL,
     "eqs_to_sage_mat_helper(equations, cols)\n"
     "--\n"
     "\n"
     "A helper to convert equation into sagemath matrix"},
    {"eqs_to_linear_system", _PyCFunction_CAST(eqs_to_linear_system),
     METH_FASTCALL,
     "eqs_to_linear_system(equations, cols)\n"
     "--\n"
     "\n"
     "Convert equations to a linear system"},
    {NULL} /* Sentinel */
};

static struct PyModuleDef _internal = {PyModuleDef_HEAD_INIT, "_internal", NULL,
                                       -1, methods};

static PyObject *gf2_matrix_new(PyTypeObject *type,
                                PyObject *args,
                                PyObject *kwds) {
	GF2MatrixObject *self;
	self = (GF2MatrixObject *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->matrix = NULL;
	}
	return (PyObject *)self;
}

static int gf2_matrix_init(GF2MatrixObject *self,
                           PyObject *args,
                           PyObject *kwds) {
	if (!PyTuple_Check(args)) {
		PyErr_SetString(PyExc_TypeError, "Arguments must be a tuple");
		return -1;
	}
	if (PyTuple_GET_SIZE(args) != 2) {
		PyErr_SetString(PyExc_TypeError, "GF2Matrix requires two argument");
		return -1;
	}
	PyObject *ints = PySequence_Fast(PyTuple_GET_ITEM(args, 0),
	                                 "First argument expected a sequence");
	PyObject *colsObj = PyTuple_GET_ITEM(args, 1);
	if (!PyLong_Check(colsObj)) {
		PyErr_SetString(PyExc_TypeError, "Second argument must be an integer");
		return -1;
	}
	Py_ssize_t cols = PyLong_AsSsize_t(colsObj);
	if (cols == -1 && PyErr_Occurred()) {
		Py_DECREF(ints);
		return -1;
	}
	if (cols <= 0) {
		Py_DECREF(ints);
		PyErr_SetString(PyExc_ValueError, "Number of columns must be positive");
		return -1;
	}
	Py_ssize_t rows = PySequence_Fast_GET_SIZE(ints);
	if (rows <= 0) {
		Py_DECREF(ints);
		PyErr_SetString(PyExc_ValueError, "Number of rows must be positive");
		return -1;
	}
	self->matrix = mzd_init(rows, cols);
	if (self->matrix == NULL) {
		Py_DECREF(ints);
		PyErr_SetString(PyExc_RuntimeError, "Failed to allocate matrix");
		return -1;
	}
	for (Py_ssize_t i = 0; i < rows; i++) {
		PyObject *row = PySequence_Fast_GET_ITEM(ints, i);
		if (!PyLong_Check(row)) {
			Py_DECREF(ints);
			PyErr_SetString(PyExc_TypeError,
			                "Row must be a sequence of integers");
			return -1;
		}
		PyLongObject *row_long = (PyLongObject *)row;
		Iter_PyLong_Bits(row_long, self->matrix->ncols,
		                 { mzd_write_bit(self->matrix, i, bitcnt, bit); }, {});
	}
	Py_DECREF(ints);
	return 0;
}

static void gf2_matrix_dealloc(GF2MatrixObject *self) {
	PyObject_GC_UnTrack(self);
	mzd_free(self->matrix);
	PyObject_GC_Del(self);
}

PyObject *gf2_matrix_richcompare(GF2MatrixObject *self,
                                 PyObject *other,
                                 int op) {
	if (op != Py_EQ) {
		Py_RETURN_NOTIMPLEMENTED;
	}
	if (self == (GF2MatrixObject *)other) {
		return PythonTrue;
	}
	if (!PyObject_TypeCheck(other, &GF2Matrix_Type)) {
		return PythonFalse;
	}
	GF2MatrixObject *other_matrix = (GF2MatrixObject *)other;
	if (self->matrix->nrows != other_matrix->matrix->nrows ||
	    self->matrix->ncols != other_matrix->matrix->ncols) {
		return PythonFalse;
	}
	int cmp = mzd_cmp(self->matrix, other_matrix->matrix);
	return cmp == 0 ? PythonTrue : PythonFalse;
}

static int gf2_matrix_get_set_item_parse_key(GF2MatrixObject *self,
                                             PyObject *key,
                                             unsigned long *iptr,
                                             unsigned long *jptr) {
	if (!PyTuple_Check(key)) {
		PyErr_SetString(PyExc_TypeError, "Arguments must be a tuple");
		return 0;
	}
	if (PyTuple_GET_SIZE(key) != 2) {
		PyErr_SetString(PyExc_TypeError, "GF2Matrix requires two argument");
		return 0;
	}
	PyObject *i = PyTuple_GET_ITEM(key, 0);
	PyObject *j = PyTuple_GET_ITEM(key, 1);
	if (!PyLong_Check(i) || !PyLong_Check(j)) {
		PyErr_SetString(PyExc_TypeError, "Index must be an integer");
		return 0;
	}
	if (PyLong_AsLong(i) < 0 || PyLong_AsLong(i) >= self->matrix->nrows) {
		PyErr_SetString(PyExc_IndexError, "Row index out of range");
		return 0;
	}
	if (PyLong_AsLong(j) < 0 || PyLong_AsLong(j) >= self->matrix->ncols) {
		PyErr_SetString(PyExc_IndexError, "Column index out of range");
		return 0;
	}
	*iptr = PyLong_AsUnsignedLong(i);
	*jptr = PyLong_AsUnsignedLong(j);
	return 1;
}

static PyObject *gf2_matrix_getitem(GF2MatrixObject *self, PyObject *key) {
	unsigned long i, j;
	if (!gf2_matrix_get_set_item_parse_key(self, key, &i, &j)) {
		return NULL;
	}
	return mzd_read_bit(self->matrix, i, j) ? PythonTrue : PythonFalse;
}

static int *gf2_matrix_setitem(GF2MatrixObject *self,
                               PyObject *key,
                               PyObject *value) {
	unsigned long i, j;
	if (!gf2_matrix_get_set_item_parse_key(self, key, &i, &j)) {
		return NULL;
	}
	if (value == NULL) {
		PyErr_SetString(PyExc_TypeError, "Deletion is not supported");
		return NULL;
	}
	if (PyObject_IsTrue(value)) {
		mzd_write_bit(self->matrix, i, j, 1);
	} else if (PyObject_Not(value)) {
		mzd_write_bit(self->matrix, i, j, 0);
	} else {
		PyErr_SetString(PyExc_TypeError, "Value must be a boolean");
		return NULL;
	}
	return NULL;
}

static PyObject *gf2_matrix_get_rows(GF2MatrixObject *self, void *closure) {
	return PyLong_FromSsize_t(self->matrix->nrows);
}

static PyObject *gf2_matrix_get_cols(GF2MatrixObject *self, void *closure) {
	return PyLong_FromSsize_t(self->matrix->ncols);
}

static inline PyObject *mzd_vector_to_pylong2(char *buf, mzd_t *v, rci_t r) {
	rci_t len = v->ncols;
	// buf[len] must be 0
	for (int c = 0; c < len; c++) {
		buf[len - 1 - c] = mzd_read_bit(v, r, c) ? '1' : '0';
	}
	return PyLong_FromString(buf, NULL, 2);
}

static PyObject *gf2_matrix_to_list(GF2MatrixObject *self) {
	PyObject *ret = PyList_New(self->matrix->nrows);
	char *str = malloc(self->matrix->ncols + 1);
	str[self->matrix->ncols] = '\0';
	for (rci_t r = 0; r < self->matrix->nrows; r++) {
		PyList_SET_ITEM(ret, r, mzd_vector_to_pylong2(str, self->matrix, r));
	}
	free(str);
	return ret;
}

static PyObject *gf2_matrix_solve_right(GF2MatrixObject *self,
                                        PyObject *const *args,
                                        Py_ssize_t nargs) {
	GF2MatrixObject *other;
	if (nargs != 1) {
		PyErr_SetString(PyExc_TypeError, "solve_right requires 1 argument");
		return NULL;
	}
	if (!PyObject_TypeCheck(args[0], &GF2Matrix_Type)) {
		PyErr_SetString(PyExc_TypeError, "Argument must be a GF2Matrix");
		return NULL;
	}
	other = (GF2MatrixObject *)args[0];
	if (other->matrix->nrows != self->matrix->nrows) {
		PyErr_SetString(PyExc_ValueError, "Number of rows must be equal");
		return NULL;
	}
	mzd_t *A = mzd_copy(NULL, self->matrix);
	mzd_t *B = mzd_copy(NULL, other->matrix);
	if (!mzd_solve_left(A, B, 0, 1)) {
		B->nrows = self->matrix->ncols;
		GF2MatrixObject *ret =
		    PyObject_GC_New(GF2MatrixObject, &GF2Matrix_Type);
		ret->matrix = mzd_transpose(NULL, B);
		PyObject_GC_Track(ret);
		mzd_free(A);
		mzd_free(B);
		return (PyObject *)ret;
	}
	mzd_free(A);
	mzd_free(B);
	return PythonNone;
}

static PyObject *gf2_matrix_right_kernel(GF2MatrixObject *self,
                                         PyObject *const *args,
                                         Py_ssize_t nargs) {
	if (nargs != 0) {
		PyErr_SetString(PyExc_TypeError, "right_kernel requires 0 argument");
		return NULL;
	}
	mzd_t *A = mzd_copy(NULL, self->matrix);
	mzd_t *ker = mzd_kernel_left_pluq(A, 0);
	mzd_free(A);
	if (ker == NULL) {
		return PythonNone;
	}
	GF2MatrixObject *ret = PyObject_GC_New(GF2MatrixObject, &GF2Matrix_Type);
	ret->matrix = mzd_transpose(NULL, ker);
	PyObject_GC_Track(ret);
	mzd_free(ker);
	return (PyObject *)ret;
}

static PyObject *gf2_matrix_multiply(GF2MatrixObject *self, PyObject *arg) {
	GF2MatrixObject *other;
	if (!PyObject_TypeCheck(arg, &GF2Matrix_Type)) {
		PyErr_SetString(PyExc_TypeError, "Argument must be a GF2Matrix");
		return NULL;
	}
	other = (GF2MatrixObject *)arg;
	if (self->matrix->ncols != other->matrix->nrows) {
		PyErr_SetString(
		    PyExc_ValueError,
		    "Number of columns of the first matrix must be equal to "
		    "the number of rows of the second matrix");
		return NULL;
	}
	mzd_t *result = mzd_mul(NULL, self->matrix, other->matrix, 0);
	GF2MatrixObject *ret = PyObject_GC_New(GF2MatrixObject, &GF2Matrix_Type);
	ret->matrix = result;
	PyObject_GC_Track(ret);
	return (PyObject *)ret;
}

static PyMappingMethods gf2_matrix_mapping = {
    .mp_subscript = (binaryfunc)gf2_matrix_getitem,
    .mp_ass_subscript = (objobjargproc)gf2_matrix_setitem,
};

static PyGetSetDef gf2_matrix_getsetdef[] = {
    {"rows", (getter)gf2_matrix_get_rows, NULL, "Number of rows", NULL},
    {"cols", (getter)gf2_matrix_get_cols, NULL, "Number of columns", NULL},
    {NULL} /* Sentinel */
};

static PyMethodDef gf2_matrix_methods[] = {
    {"to_list", _PyCFunction_CAST(gf2_matrix_to_list), METH_NOARGS,
     "to_list(s)\n"
     "--\n"
     "\n"
     "Convert matrix to list"},
    {"solve_right", _PyCFunction_CAST(gf2_matrix_solve_right), METH_FASTCALL,
     "solve_right(other)\n"
     "--\n"
     "\n"
     "Solve a linear system"},
    {"right_kernel", _PyCFunction_CAST(gf2_matrix_right_kernel), METH_FASTCALL,
     "right_kernel()\n"
     "--\n"
     "\n"
     "Compute right kernel"},
    {NULL} /* Sentinel */
};

static PyNumberMethods gf2_matrix_numbers = {
    .nb_matrix_multiply = (binaryfunc)gf2_matrix_multiply,
};

static PyTypeObject GF2Matrix_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name = "_internal.GF2Matrix",
    .tp_basicsize = sizeof(AffineSpaceIteratorObject),
    .tp_itemsize = 0,
    .tp_new = gf2_matrix_new,
    .tp_init = (initproc)gf2_matrix_init,
    .tp_dealloc = (destructor)gf2_matrix_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_traverse = (traverseproc)nop,
    .tp_clear = (inquiry)nop,
    .tp_doc = NULL,
    .tp_as_mapping = &gf2_matrix_mapping,
    .tp_getset = gf2_matrix_getsetdef,
    .tp_methods = gf2_matrix_methods,
    .tp_as_number = &gf2_matrix_numbers,
    .tp_richcompare = (richcmpfunc)gf2_matrix_richcompare,
};

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
	INIT_TYPE(AffineSpaceIterator_Type);
	INIT_TYPE(AffineSpaceIteratorSlow_Type);
	INIT_TYPE(GF2Matrix_Type);
	PyObject *mod = PyModule_Create(&_internal);
	if (mod == NULL)
		return NULL;
	ADD_TYPE(mod, AffineSpace_Type);
	ADD_TYPE(mod, AffineSpaceIterator_Type);
	ADD_TYPE(mod, AffineSpaceIteratorSlow_Type);
	ADD_TYPE(mod, GF2Matrix_Type);
	return mod;
}
