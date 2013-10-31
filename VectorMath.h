
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

static int VectorObject_cmp(VectorObject *self, VectorObject *other); // declared differently because it is set as a cmp function instead of a class method
static PyObject *VectorObject_copy(VectorObject *self);
static PyObject *VectorObject_add(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_sub(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_length(VectorObject *self);
static PyObject *VectorObject_crossProduct(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_dotProduct(VectorObject *self, VectorObject *other);
static PyObject *VectorObject_angle(VectorObject *self, VectorObject *other); 

static char module_docstring[] =
    "This module provides an interface for vector processing in C.";
static char copy_docstring[] =
	"Returns a copy of the vector argument.";
static char add_docstring[] =
	"Adds vector A and vector B.";
static char sub_docstring[] =
	"Subtracts vector B from vector A.";
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
