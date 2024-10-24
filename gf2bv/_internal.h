#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <m4ri/m4ri.h>

typedef struct {
	PyObject_HEAD mzd_t *origin;
	mzd_t *basis;
} AffineSpaceObject;

typedef struct {
	PyObject_HEAD AffineSpaceObject *space;
	// uint8_t *state;
	union {
		uint8_t *state;  // for naive enumeration
		struct {
			uint64_t idx;
			mzd_t *cur;
		} gray;  // for gray code enumeration
	};
	char *str;
} AffineSpaceIteratorObject;

typedef struct {
	PyObject_HEAD mzd_t *matrix;
} GF2MatrixObject;

#define SOLVE_MODE_SINGLE 0
#define SOLVE_MODE_AFFINE_SPACE 1

PyObject *m4ri_solve(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *to_bits(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *mul_bit_quad(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *xor_list(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *list_where(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *eqs_to_sage_mat_helper(PyObject *self,
                                 PyObject *const *args,
                                 Py_ssize_t nargs);
