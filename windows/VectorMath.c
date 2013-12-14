# include <Python.h>
# include <structmember.h>
# include "VectorMath.h"
# include "backend\\vector.h"

/**
 * Following code is defining the Vector object that VectorMath operates on.
 * 
 * It is constructed using C ints and C doubles, and has getters and setters
 * that allow interfacing with them in order to avoid constant conversion
 * back and forth everytime an operation is performed. This way we only have to 
 * convert the C values to Python objects, and only do it when we need to access
 * the members from Python
 */
 
/*
 ***********************************************************************
 * 
 * Define methods of Vector object
 * 
 ***********************************************************************
 */

static PyMethodDef PyVector_methods[] = {
    {"copy", PyVector_copy, METH_NOARGS, copy_docstring},
    {"length", PyVector_length, METH_NOARGS, magnitude_docstring},
    {"magnitude", PyVector_length, METH_NOARGS, magnitude_docstring},
    {"cross", PyVector_crossProduct, METH_O, crossProduct_docstring},
    {"crossProduct", PyVector_crossProduct, METH_O, crossProduct_docstring},
    {"dot", PyVector_dotProduct, METH_O, dotProduct_docstring},
    {"dotProduct", PyVector_dotProduct, METH_O, dotProduct_docstring},
    {"angle", PyVector_angle, METH_O, angle_docstring}, 
    {"angularDifference", PyVector_angle, METH_O, angle_docstring}, 
    {"norm", PyVector_normalize, METH_NOARGS, NULL},
    {"normalize", PyVector_normalize, METH_NOARGS, NULL},
    {"rotate", PyVector_rotate, METH_O, NULL},
    {NULL, NULL, 0, NULL}
}; 

static PyNumberMethods PyVector_as_number = {
    (binaryfunc)PyVector_add,        /*nb_add*/
    (binaryfunc)PyVector_sub,        /*nb_subtract*/
    (binaryfunc)PyVector_mul,        /*nb_multiply*/
    (binaryfunc)PyVector_div,        /*nb_divide*/
    0,                       /*nb_remainder*/
    0,                       /*nb_divmod*/
    0,                       /*nb_power*/
    (unaryfunc)PyVector_neg,        /*nb_negative*/
    0,                       /*nb_positive*/
    (unaryfunc)PyVector_length,        /*nb_absolute*/
    (inquiry)PyVector_nonzero,       /*nb_nonzero*/
    (unaryfunc)PyVector_neg,           /*nb_invert*/
    0,                          /*nb_lshift*/
    0,                          /*nb_rshift*/
    0,                          /*nb_and*/
    0,                          /*nb_xor*/
    0,                          /*nb_or*/
    0,                          /*nb_coerce*/
    0,                          /*nb_int*/
    0,                          /*nb_long*/
    0,                          /*nb_float*/
    0,                          /*nb_oct*/
    0,                          /*nb_hex*/
    0,                          /*nb_inplace_add*/
    0,                          /*nb_inplace_subtract*/
    0,                          /*nb_inplace_multiply*/
    0,                          /*nb_inplace_divide*/
    0,                          /*nb_inplace_remainder*/
    0,                          /*nb_inplace_power*/
};

static PyGetSetDef PyVector_getset[] = {
	{"components", PyVector_getComponents, PyVector_setComponents, NULL, NULL},
	{"demensions", PyVector_getDemensions, PyVector_setDemensions, NULL, NULL},
	{NULL, NULL, NULL, NULL, NULL}
};
 

static PyMethodDef module_methods[] = {
	{NULL, NULL, 0, NULL}
};
 
/*
 ***********************************************************************
 * 
 * Package everything into one object
 * 
 ***********************************************************************
 */

static PyTypeObject PyVectorType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "VectorMath.Vector",       /*tp_name*/
    sizeof(PyVector),      /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PyVector_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    PyVector_cmp,          /*tp_compare*/
    0,                         /*tp_repr*/
    &PyVector_as_number,    /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Vector quantity",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    PyVector_methods,      /* tp_methods */
    0,                         /* tp_members */
    PyVector_getset,       /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PyVector_init, /* tp_init */
    0,                           /* tp_alloc */
    PyVector_new            /* tp_new */
};

/*
 ***********************************************************************
 * 
 * Define getters and setters for members, as they are part of the 
 * C struct, and cannot be accessed through Python otherwise
 * 
 ***********************************************************************
 */
 
static PyObject *PyVector_getComponents(PyVector *self) {
	
	PyTupleObject *retList;     // Tuple containing components (faster than lists)
	int i;
	
	retList = PyTuple_New(self->v->demensions); // create new tuple
	
	if (!retList) {             // Ensure the list was created properly
		return (PyObject *) NULL;
	}
	
	for (i=0; i<self->v->demensions; i++) {   // Iterate through components and place them in the tuple
		PyTuple_SET_ITEM((PyObject *) retList, i, PyFloat_FromDouble(self->v->components[i]));
	}
	
	return (PyObject *) retList; 
}

static int PyVector_setComponents(PyVector *self, PyObject *l) {
	
	int i;
	
	if (PyTuple_Check(l)) {
		
		for (i=0; i<self->v->demensions; i++) {  // Iterate through tuple and assign components
			self->v->components[i] = PyFloat_AsDouble(PyTuple_GET_ITEM((PyObject *) l, i));
		}
		
	} else if (PyList_Check(l)) {

		for (i=0; i<self->v->demensions; i++) {  // Iterate through list and assign components
			self->v->components[i] = PyFloat_AsDouble(PyList_GET_ITEM((PyObject *) l, i));
		}
		
	} else {
		return -1;
	}
	
	return 0;
}

static PyObject *PyVector_getDemensions(PyVector *self) {
	
	PyLongObject *n;     // Tuple containing components (faster than lists)

	n = PyLong_FromLong((long) self->v->demensions);
	
	return (PyObject *) n; 
}

static int PyVector_setDemensions(PyVector *self, PyObject *i) {
	
	self->v->demensions = (int) PyLong_AsLong(i);  // attempt to read int
	
	return 0;
}
 
/*
 ***********************************************************************
 * 
 * Define constructors and destructors
 *  - dealloc
 *  - new
 *  - init
 * 
 ***********************************************************************
 */

static PyObject *PyVector_new(PyTypeObject *type, PyObject *args) {
	
    PyObject *newObject;                 // Vector object being created
    
    if (type != &PyVectorType) { 
		assert(PyType_IsSubtype(type, &PyVectorType));  // Make sure it's a subclass of a vector
	}
	
    newObject = type->tp_alloc(type, 0);                // allocate type
    
    if (!newObject) {      
		return (PyObject *)NULL;
	}
	
	((PyVector *)newObject)->v = malloc(sizeof(Vector));  // allocate C-vector member

    return newObject;
}

static int PyVector_init(PyVector *self, PyObject *args) {
	
	PyObject *compList;
	PyObject *(*getitem)(PyObject *, int);
	int i;

    if (!PyArg_ParseTuple(args, "O", &compList)) { // Extract demensions tuple
        return -1; 
	}
	
	if (PyTuple_Check(compList)) {
		
		getitem = PyTuple_GetItem;
		self->v->demensions = PyTuple_GET_SIZE(compList);      // Get number of demensions
		
	} else if (PyList_Check(compList)) {
		
		getitem = PyList_GetItem;
		self->v->demensions = PyList_GET_SIZE(compList);      // Get number of demensions
		
	} else {
		
		return -1;
		
	}
	
	if (self->v->demensions > 3) {          // cap the number of components at 3. Not sure if it should spit back an error
		self->v->demensions = 3;
	}
	
	for (i=0; i<self->v->demensions; i++) {  // Iterate through tuple and set components
		self->v->components[i] = PyFloat_AsDouble(getitem(compList, i));
	}
	
    return 0;       // Return success
}

static void PyVector_dealloc(PyVector *self) {
	
	free(self->v);   // Free initalized vector
    self->ob_type->tp_free((PyObject*)self);  // Free PyObject
    
}

/*
 ***********************************************************************
 * 
 * Core functions of VectorMath library
 * 
 ***********************************************************************
 */
 
static PyObject *PyVector_copy(PyVector *self) {
	
	PyVector *result = (PyVector *) PyVector_new(self->ob_type, NULL);
	Vector *temp, *v = vectorCopy(self->v);
	
	temp = result->v;
	result->v = v;
	
	free(temp);
	
	return (PyObject *) result;
}
 
static PyObject *PyVector_add(PyVector *self, PyVector *other) {
	
	PyVector *result;    // Object being returned
	Vector *vr, *temp;  // C-vectors of a, b, and return value
	
	vr = vectorAdd(self->v, other->v);  // Add the two vectors and get the new one in return
	 
	if (!vr) {    // Ensure we recieved a vector 
		return NULL;
	}
	
	result    = (PyVector *) PyVector_new(self->ob_type, NULL); // Create new PyVector
	temp      = result->v;      // this whole process basically frees the vector currently allocated, and sets the pointer equal to the resultant vector
	result->v = vr;	            // It's faster than copying the data over because either way we have to free a vector. So this way we just toss a pointer around instead of 3 doubles
	
	vectorFree(temp);
	
	return (PyObject *) result;
}

static PyObject *PyVector_sub(PyVector *self, PyVector *other) {
	
	PyVector *result;    // Object being returned
	Vector *vr, *temp;  // C-vectors of a, b, and return value
	
	vr = vectorSub(self->v, other->v);  // Add the two vectors and get the new one in return
	 
	if (!vr) {    // Ensure we recieved a vector 
		return NULL;
	}
	
	result    = (PyVector *) PyVector_new(self->ob_type, NULL); // Create new PyVector
	temp      = result->v;
	result->v = vr;	
	
	vectorFree(temp);
	
	return (PyObject *) result;
}

static PyObject *PyVector_crossProduct(PyVector *self, PyVector *other) {
	
	PyVector *retObject;
	Vector *vr, *temp;
	
	vr = crossProduct(self->v, other->v); // Calculate cross-product
	
	if (!vr) { // ensure we have a vector (High probability that we don't since cross product only works on 3D vectors)
		return NULL;
	}
	
	retObject = (PyVector *) PyVector_new(self->ob_type, NULL);
	temp      = retObject->v;
	retObject->v = vr;
	
	vectorFree(temp);
	
	return (PyObject *) retObject;
}

static PyObject *PyVector_mul(PyVector *self, PyObject *other) {
	
	double k;
	PyVector *ret;
	Vector *temp, *vr;
	
	if (PyObject_TypeCheck(other, self->ob_type)) {  // If the other object is a vector
	
		return PyVector_dotProduct(self, other);     // make the crude assumption they meant dot product.
		 
	} else {
		
		k   = PyFloat_AsDouble(other);          // Otherwise, we assume it's a number
		ret = (PyVector *) PyVector_new(self->ob_type, NULL);  // Create new pyvector object
		
		vr     = vectorMul(self->v, k);       // multiply the vector by constant k
		temp   = ret->v;              // set temp equal to currently allocated vector
		ret->v = vr;                  // set the vector pointer to the resulting vector
		
		vectorFree(temp);             // free the previously allocated vector
		
	    return (PyObject *) ret;
	}
}

static PyObject *PyVector_div(PyVector *self, PyObject *other) {
	
	double k;
	PyVector *ret;
	Vector *temp, *vr;
	
	if (PyObject_TypeCheck(other, self->ob_type)) {  // If the other object is a vector
	
		return (PyObject *)NULL; // have no defined way to divide vectors. 
		 
	} else {  // We can however divide by scalars!
		
		k   = PyFloat_AsDouble(other);          // Otherwise, we assume it's either a float or an int
		ret = (PyVector *) PyVector_new(self->ob_type, NULL);  // Create new pyvector object
		
		vr     = vectorDiv(self->v, k);       // multiply the vector by constant k
		temp   = ret->v;              // set temp equal to currently allocated vector
		ret->v = vr;                  // set the vector pointer to the resulting vector
		
		vectorFree(temp);             // free the previously allocated vector
		
	    return (PyObject *) ret;
	}
}

static PyObject *PyVector_neg(PyVector *self) {
	
	PyVector *ret;
	
	ret = PyVector_new(self->ob_type, NULL); // create PyVector object
	
	free(ret->v);
	
	ret->v = vectorCopy(self->v);            // copy c-vector member
	
	ret->v->components[0] = -ret->v->components[0];  // negate each component
	ret->v->components[1] = -ret->v->components[1];
	ret->v->components[2] = -ret->v->components[2];
	
	return (PyObject *) ret;
}

static PyObject *PyVector_normalize(PyVector *self) {
	
	PyVector *newVect;
	
	newVect = (PyVector *) PyVector_new(self->ob_type, NULL);
	
	free(newVect->v);
	
	newVect->v = vectorNormalize(self->v);

	return (PyObject *) newVect;
}

static PyObject *PyVector_length(PyVector* self) {
	
    PyObject *result;
    double length;
   
    length = vectorLength(self->v);
    
    result = PyFloat_FromDouble(length);
    
    return result;
}

static PyObject *PyVector_dotProduct(PyVector *self, PyVector *other) {
	
	PyObject *retObject;
	double result;
	
	result = dotProduct(self->v, other->v); // Calculate scalar-product
	
	retObject = PyFloat_FromDouble(result);
	
	return retObject;
}

static PyObject *PyVector_angle(PyVector *self, PyVector *other) {
	
	PyObject *retObject;
	double result;
	
	result = angle(self->v, other->v); // Calculate scalar-product
	
	retObject = PyFloat_FromDouble(result);
	
	return retObject;
}

static PyObject *PyVector_rotate(PyVector *self, PyVector *rotation) {
	
	PyVector *newVect;
	
	newVect = (PyVector *) PyVector_new(self->ob_type, NULL);
	
	free(newVect->v);
	
	newVect->v = vectorRotate(self->v, rotation->v);

	return (PyObject *) newVect;
}

static int PyVector_cmp(PyVector *self, PyVector *other) {
	
	if (!PyObject_TypeCheck(other, &PyVectorType)) { 
		// Ensure the object is a vector
		return -2;
	} else if (self->v->demensions < other->v->demensions) {
		// If vector A has less demensions than vector B, we're gonna say A < B
		return -1;
	} else if (self->v->demensions > other->v->demensions) {
		// If B has less demensions, then we say A > B
		return 1;
	}  else {
		// If we have the same number of demensions, then we compare magnitudes
		return vectorCompare(self->v, other->v);
	}
}

static int PyVector_nonzero(PyVector *self) {
	
	if (vectorLength(self->v)) { // If it has any sort of length, then it is not zero
		return 1;
	}
	
	return 0;  // Otherwise it is zero
}

/*
 ***********************************************************************
 * 
 * Initialize module
 * 
 ***********************************************************************
 */

PyMODINIT_FUNC initVectorMath(void) {
    PyObject *m = Py_InitModule3("VectorMath", module_methods, module_docstring); // Setup module
    if (m == NULL)
        return;
        
    if (PyType_Ready(&PyVectorType) < 0) {     // Finalize type
		return;
	}
        
    Py_INCREF(&PyVectorType);                  // Increment reference to pass to Python
    PyModule_AddObject(m, "Vector", (PyObject *)&PyVectorType);  // Add to objects list so we can see it
}
