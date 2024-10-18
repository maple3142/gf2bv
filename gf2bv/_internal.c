#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <m4ri/m4ri.h>

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

static PyObject *mzd_vector_to_pylong(char *buf, mzd_t *v) {
	rci_t len = v->ncols;
	// buf[len] must be 0
	for (int c = 0; c < len; c++) {
		buf[len - 1 - c] = mzd_read_bit(v, 0, c) ? '1' : '0';
	}
	return PyLong_FromString(buf, NULL, 2);
}

// iterator reference: https://github.com/Lydxn/xorsat/blob/789fed013292f060c026be8d1990631041969e40/xorsat/_xorsatmodule.c#L1377
typedef struct {
	PyObject_HEAD mzd_t *sol0;
	mzd_t *kernel;
	uint8_t *state;
	char *str;
} SolutionIterObject;

static PyObject *solveiter_next(SolutionIterObject *it) {
	// TODO: use gray code for more efficient enumeration

	mzd_t *result;
	rci_t n, r;
	int sentinel;

	n = it->kernel->nrows;
	if (it->state[n])
		return NULL; /* StopIteration */

	result = mzd_copy(NULL, it->sol0);
	for (r = 0; r < n; r++)
		if (it->state[r])
			mzd_combine_even_in_place(result, 0, 0, it->kernel, r, 0);

	sentinel = 1;
	for (r = 0; r < n; r++) {
		if ((it->state[r] ^= 1)) {
			sentinel = 0;
			break;
		}
	}
	it->state[n] = sentinel;

	PyObject *ret = mzd_vector_to_pylong(it->str, result);
	mzd_free(result);
	return ret;
}

static void solveiter_dealloc(SolutionIterObject *self) {
	PyObject_GC_UnTrack(self);
	mzd_free(self->sol0);
	mzd_free(self->kernel);
	free(self->state);
	free(self->str);
	PyObject_GC_Del(self);
}

static PyTypeObject SolutionIter_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name =
        "_internal.SolutionIter",
    .tp_basicsize = sizeof(SolutionIterObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)solveiter_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_doc = NULL,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)solveiter_next,
};

static PyObject *m4ri_solve(PyObject *self,
                            PyObject *const *args,
                            Py_ssize_t nargs) {
	PyObject *linsys_list;
	Py_ssize_t cols;
	int all = 0;
	// parse the arguments: (list, cols, all)
	// if (!PyArg_ParseTuple(args, "O!np", &PyList_Type, &linsys_list, &cols,
	//                       &all))
	// 	return NULL;
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
	if (cols < 0) {
		if (cols == -1 && PyErr_Occurred()) {
			return NULL;
		}
		PyErr_SetString(PyExc_ValueError,
		                "Number of columns must be non-negative");
		return NULL;
	}
	all = PyObject_IsTrue(args[2]);
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
		Py_ssize_t n_digits = PyLong_DigitCount(v);
		Py_ssize_t c = -1;
		for (Py_ssize_t i = 0; i < n_digits && c < cols; i++) {
			digit d = GET_OB_DIGITS(v)[i];
			for (int j = 0; j < PyLong_SHIFT && c < cols; j++) {
				if (c == -1) {
					mzd_write_bit(B, r, 0, d & 1);
					B_is_not_zero |= d & 1;
				} else {
					mzd_write_bit(A, r, c, d & 1);
				}
				d >>= 1;
				c++;
			}
		}
	}

	// first, find the base solution
	mzd_t *sol0;
	{
		// for kernel, only needed if all solutions are needed and B != 0
		mzd_t *A_copy = 0;
		if (B_is_not_zero) {
			if (all) {
				A_copy = mzd_copy(NULL, A);
			}
			// printf("A->nrows: %d, A->ncols: %d\n", A->nrows, A->ncols);
			// printf("B->nrows: %d, B->ncols: %d\n", B->nrows, B->ncols);
			// this is actually solve_right in Sage...
			if (mzd_solve_left(A, B, 0, 1) != 0) {
				mzd_free(A);
				mzd_free(B);
				if (all) {
					mzd_free(A_copy);
				}
				return Py_None;
			}
			// the base solution is stored in B as column vector
			B->nrows = cols;
			sol0 = mzd_transpose(0, B);
			mzd_free(A);  // A is overwritten by mzd_solve_left
			A = A_copy;
		} else {
			// the trivial solution
			sol0 = mzd_init(1, cols);
		}
	}
	mzd_free(B);
	// if all: A, sol0 is valid
	// else: A is valid or 0, sol0 is valid

	// if we only need one solution, return it
	if (!all) {
		char *str = malloc(cols + 1);
		str[cols] = '\0';
		PyObject *ret = mzd_vector_to_pylong(str, sol0);
		free(str);
		mzd_free(sol0);
		if (A) {
			mzd_free(A);
		}
		return ret;
	}

	// if all solutions are needed, we need to find the kernel

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
	// now, only sol0 and tker is valid

	SolutionIterObject *it =
	    PyObject_GC_New(SolutionIterObject, &SolutionIter_Type);
	it->sol0 = sol0;
	it->kernel = tker;
	it->state = malloc(tker->nrows + 1);
	it->str = malloc(cols + 1);
	it->str[cols] = '\0';
	memset(it->state, 0, tker->nrows + 1);

	return (PyObject *)it;
}

#if PY_VERSION_HEX >= 0x030C0000
// 3.12 and later they are immortal: https://peps.python.org/pep-0683/
#define PythonTrue Py_True
#define PythonFalse Py_False
#else
// pre 3.12
#define PythonTrue Py_NewRef(Py_True)
#define PythonFalse Py_NewRef(Py_False)
#endif

static PyObject *to_bits(PyObject *self,
                         PyObject *const *args,
                         Py_ssize_t nargs) {
	// parse the arguments: (n: int, a: bigint)
	// convert a bigint to a list of booleans of length n, little-endian
	Py_ssize_t n;
	PyObject *a;
	// if (!PyArg_ParseTuple(args, "nO!", &n, &PyLong_Type, &a))
	// 	return NULL;
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
	Py_ssize_t n_digits_a = PyLong_DigitCount(a_long);
	Py_ssize_t c = 0;
	for (Py_ssize_t i = 0; i < n_digits_a && c < n; i++) {
		digit d = GET_OB_DIGITS(a_long)[i];
		for (int j = 0; j < PyLong_SHIFT && c < n; j++) {
			PyList_SetItem(list, c, d & 1 ? PythonTrue : PythonFalse);
			d >>= 1;
			c++;
		}
	}
	for (; c < n; c++) {
		PyList_SetItem(list, c, PythonFalse);
	}
	return list;
}

void set_bits(char *bits, Py_ssize_t n, PyLongObject *o) {
	Py_ssize_t n_digits_a = PyLong_DigitCount(o);
	Py_ssize_t c = 0;
	for (Py_ssize_t i = 0; i < n_digits_a && c < n; i++) {
		digit d = GET_OB_DIGITS(o)[i];
		for (int j = 0; j < PyLong_SHIFT && c < n; j++) {
			bits[c] = d & 1 ? 1 : 0;
			d >>= 1;
			c++;
		}
	}
	for (; c < n; c++) {
		bits[c] = 0;
	}
}

static PyObject *mul_bit_quad(PyObject *self,
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
     "Solve a linear system over GF(2) with M4RI"},
    {"to_bits", _PyCFunction_CAST(to_bits), METH_FASTCALL,
     "Convert an integer to a list of bits"},
    {"mul_bit_quad", _PyCFunction_CAST(mul_bit_quad), METH_FASTCALL,
     "Multiply two linear symbolic bits to a linearized quadratic symbolic "
     "bit"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef _internal = {PyModuleDef_HEAD_INIT, "_internal", NULL,
                                       -1, methods};

PyMODINIT_FUNC PyInit__internal(void) {
	return PyModule_Create(&_internal);
}
