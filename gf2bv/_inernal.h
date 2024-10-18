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
	uint8_t *state;
	char *str;
} AffineSpaceIterableObject;

#define SOLVE_MODE_SINGLE 0
#define SOLVE_MODE_ALL 1

PyObject *m4ri_solve(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *to_bits(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
PyObject *mul_bit_quad(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
