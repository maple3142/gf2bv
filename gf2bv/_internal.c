#include "_internal.h"

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
	INIT_TYPE(AffineSpaceIterator_Type);
	INIT_TYPE(AffineSpaceIteratorSlow_Type);
	PyObject *mod = PyModule_Create(&_internal);
	if (mod == NULL)
		return NULL;
	ADD_TYPE(mod, AffineSpace_Type);
	ADD_TYPE(mod, AffineSpaceIterator_Type);
	ADD_TYPE(mod, AffineSpaceIteratorSlow_Type);
	return mod;
}
