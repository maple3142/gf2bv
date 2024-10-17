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

// iterator reference: https://github.com/Lydxn/xorsat/blob/789fed013292f060c026be8d1990631041969e40/xorsat/_xorsatmodule.c#L1377
typedef struct {
	PyObject_HEAD rci_t cols;
	mzd_t *sol0;
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

	for (rci_t c = 0; c < it->cols; c++) {
		it->str[it->cols - 1 - c] = mzd_read_bit(result, 0, c) ? '1' : '0';
	}
	mzd_free(result);
	PyObject *ret = PyLong_FromString(it->str, NULL, 2);
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
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name = "pym4ri.SolutionIter",
    .tp_basicsize = sizeof(SolutionIterObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)solveiter_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .tp_doc = NULL,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)solveiter_next,
};

static PyObject *solve(PyObject *self, PyObject *args) {
	PyObject *linsys_list;
	Py_ssize_t cols;
	int all = 0;
	// parse the arguments: (list, cols, all)
	if (!PyArg_ParseTuple(args, "O!np", &PyList_Type, &linsys_list, &cols,
	                      &all))
		return NULL;
	if (cols < 0) {
		PyErr_SetString(PyExc_ValueError,
		                "Number of columns must be non-negative");
		return NULL;
	}
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

	mzd_t *A_copy = 0;  // for kernel, only needed if all solutions are needed
	if (all) {
		A_copy = mzd_copy(NULL, A);
	}

	mzd_t *sol0;
	if (B_is_not_zero) {
		// this is actually solve_right in Sage...
		if (mzd_solve_left(A, B, 0, 1) != 0) {
			mzd_free(A);
			mzd_free(B);
			return Py_None;
		}
		// the base solution, stored in B as column vector
		sol0 = mzd_transpose(0, B);
		mzd_free(A);
		mzd_free(B);
	} else {
		mzd_free(A);
		mzd_free(B);
		sol0 = mzd_init(1, cols);  // the trivial solution
	}

	if (!all) {
		char *str = malloc(cols + 1);
		for (int c = 0; c < cols; c++) {
			str[cols - 1 - c] = mzd_read_bit(sol0, 0, c) ? '1' : '0';
		}
		str[cols] = '\0';

		PyObject *ret = PyLong_FromString(str, NULL, 2);
		free(str);
		mzd_free(sol0);
		return ret;
	}

	// and this is right_kernel_matrix
	mzd_t *ker = mzd_kernel_left_pluq(A_copy, 0);
	mzd_free(A_copy);
	mzd_t *tker;
	if (ker == NULL) {
		tker = mzd_init(0, 0);
	} else {
		tker = mzd_transpose(0, ker);
		mzd_free(ker);
	}

	SolutionIterObject *it =
	    PyObject_GC_New(SolutionIterObject, &SolutionIter_Type);
	it->cols = cols;
	it->sol0 = sol0;
	it->kernel = tker;
	it->state = malloc(tker->nrows + 1);
	it->str = malloc(cols + 1);
	it->str[cols] = '\0';
	memset(it->state, 0, tker->nrows + 1);

	return (PyObject *)it;
}

static PyMethodDef methods[] = {{"solve", solve, METH_VARARGS,
                                 "Solve a linear system over GF(2) with M4RI"},
                                {NULL, NULL, 0, NULL}};

static struct PyModuleDef pym4ri = {PyModuleDef_HEAD_INIT, "pym4ri", NULL, -1,
                                    methods};

PyMODINIT_FUNC PyInit_pym4ri(void) {
	return PyModule_Create(&pym4ri);
}
