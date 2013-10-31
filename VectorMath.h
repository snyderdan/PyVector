
# ifndef PYVECTOR_H

# define PYVECTOR_H

# ifdef __cplusplus

extern "C" {
	
# endif

# include "backend//vector.h"

typedef struct {
    PyObject_HEAD       // Python object headers
    Vector *v;
} VectorObject;

// comparisons
static int VectorObject_cmp(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_nonzero(VectorObject *self);
// vector return values
static PyObject *VectorObject_add(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_sub(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_mul(VectorObject *self, PyObject *other);
static PyObject *VectorObject_div(VectorObject *self, PyObject *other);
static PyObject *VectorObject_crossProduct(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_normalize(VectorObject *self);
static PyObject *VectorObject_copy(VectorObject *self);
// scalar return values
static PyObject *VectorObject_dotProduct(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_angle(VectorObject *self, VectorObject *other); 
static PyObject *VectorObject_neg(VectorObject *self);
static PyObject *VectorObject_length(VectorObject *self);

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
