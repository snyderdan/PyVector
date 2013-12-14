
# ifndef PYVECTOR_H

# define PYVECTOR_H

# ifdef __cplusplus

extern "C" {
	
# endif

# include "backend/vector.h"

typedef struct {
    PyObject_HEAD       // Python object headers
    Vector *v;
} PyVector;

static PyObject *PyVector_new(PyTypeObject *type, PyObject *args);
static int PyVector_init(PyVector *self, PyObject *args);
static void PyVector_dealloc(PyVector *self);
static int PyVector_setDemensions(PyVector *self, PyObject *i);
static int PyVector_setComponents(PyVector *self, PyObject *i);
static PyObject *PyVector_getComponents(PyVector *self);
static PyObject *PyVector_getDemensions(PyVector *self);

// comparisons
static int PyVector_cmp(PyVector *self, PyVector *other);
static int PyVector_nonzero(PyVector *self);
// operations returning vectors
static PyObject *PyVector_add(PyVector *self, PyVector *other);
static PyObject *PyVector_sub(PyVector *self, PyVector *other);
static PyObject *PyVector_mul(PyVector *self, PyObject *other);
static PyObject *PyVector_div(PyVector *self, PyObject *other);
static PyObject *PyVector_crossProduct(PyVector *self, PyVector *other);
static PyObject *PyVector_rotate(PyVector *self, PyVector *rotation);
static PyObject *PyVector_normalize(PyVector *self);
static PyObject *PyVector_neg(PyVector *self);
static PyObject *PyVector_copy(PyVector *self);
// operations returning scalars
static PyObject *PyVector_dotProduct(PyVector *self, PyVector *other);
static PyObject *PyVector_angle(PyVector *self, PyVector *other);
static PyObject *PyVector_length(PyVector *self);

static char module_docstring[] =
    "This module provides an interface for vector processing in C.";
static char copy_docstring[] =
	"Returns a copy of the vector argument.";
static char magnitude_docstring[] = 
	"Computes the magnitude of a vector.";
static char crossProduct_docstring[] =
	"Computes the cross-product of two vectors.";
static char dotProduct_docstring[] =
	"Computes the dot-product of two vectors.";
static char angle_docstring[] = 
	"Returns the angle between two vector quantities.";

#ifdef __cplusplus

}
#endif

#endif // PYVECTOR_H
