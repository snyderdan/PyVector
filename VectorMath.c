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
 * Define getters and setters for members, as they are part of the 
 * C struct, and cannot be accessed through Python otherwise
 * 
 ***********************************************************************
 */
 
static PyObject *VectorObject_getComponents(VectorObject *self) {
	
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

static int VectorObject_setComponents(VectorObject *self, PyObject *l) {
	
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

static PyObject *VectorObject_getDemensions(VectorObject *self) {
	
	PyIntObject *n;     // Tuple containing components (faster than lists)

	n = PyInt_FromLong((long) self->v->demensions);
	
	return (PyObject *) n; 
}

static void VectorObject_setDemensions(VectorObject *self, PyIntObject *i) {
	
	self->v->demensions = (int) PyInt_AsLong(i);
	
}
 
static PyGetSetDef VectorObject_getset[] = {
	{"components", VectorObject_getComponents, VectorObject_setComponents, NULL, NULL},
	{"demensions", VectorObject_getDemensions, VectorObject_setDemensions, NULL, NULL},
	{NULL, NULL, NULL, NULL, NULL}
};
 
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

static void VectorObject_dealloc(VectorObject *self) {
	
	vectorFree(self->v);   // Free initalized vector
    self->ob_type->tp_free((PyObject*)self);  // Free PyObject
    
}

static PyObject * VectorObject_new(PyTypeObject *type, PyObject *args) {
	
    VectorObject *self;                 // Vector object being created

    self = (VectorObject *)type->tp_alloc(type, 0); // Allocate size of vector
    
    if (!self) {         // Ensure memory was allocated
		return (PyObject *)NULL;
	}
	
	self->v = calloc(sizeof(Vector), 0);  // Allocate vector without intializing components 

    return (PyObject *)self;
}

static int VectorObject_init(VectorObject *self, PyObject *args) {
	
	PyObject *compList;
	int i;

    if (!PyArg_ParseTuple(args, "O", &compList)) { // Extract demensions tuple
        return -1; 
	}
	
	if (PyTuple_Check(compList)) {   // If we recieved tuple from Python
		
		self->v->demensions = PyTuple_GET_SIZE(compList);      // Get number of demensions
		
		for (i=0; i<self->v->demensions; i++) {  // Iterate through tuple and set components
			self->v->components[i] = PyFloat_AsDouble(PyTuple_GET_ITEM(compList, i));
		}
		
	} else if (PyList_Check(compList)) { // If we recieved a list from Python
		
		self->v->demensions = PyList_GET_SIZE(compList);      // Get number of demensions

		for (i=0; i<self->v->demensions; i++) {  // Iterate through list and set components
			self->v->components[i] = PyFloat_AsDouble(PyList_GET_ITEM(compList, i));
		}
		
	} else {
		
		return -1;
		
	}
	
    return 0;       // Return success
}

/*
 ***********************************************************************
 * 
 * Define methods of Vector object
 * 
 ***********************************************************************
 */

static PyMethodDef VectorObject_methods[] = {
    {"copy", VectorObject_copy, METH_NOARGS, copy_docstring},
    {"length", VectorObject_length, METH_NOARGS, magnitude_docstring},
    {"magnitude", VectorObject_length, METH_NOARGS, magnitude_docstring},
    {"cross", VectorObject_crossProduct, METH_O, crossProduct_docstring},
    {"crossProduct", VectorObject_crossProduct, METH_O, crossProduct_docstring},
    {"dot", VectorObject_dotProduct, METH_O, dotProduct_docstring},
    {"dotProduct", VectorObject_dotProduct, METH_O, dotProduct_docstring},
    {"angle", VectorObject_angle, METH_O, angle_docstring}, 
    {"angularDifference", VectorObject_angle, METH_O, angle_docstring}, 
    {"normalize", VectorObject_normalize, METH_O, NULL},
    {NULL, NULL, 0, NULL}
}; 

static PyNumberMethods VectorObject_as_number = {
    (binaryfunc)VectorObject_add,        /*nb_add*/
    (binaryfunc)VectorObject_sub,        /*nb_subtract*/
    (binaryfunc)VectorObject_mul,        /*nb_multiply*/
    (binaryfunc)VectorObject_div,        /*nb_divide*/
    0,                       /*nb_remainder*/
    0,                       /*nb_divmod*/
    0,                       /*nb_power*/
    (unaryfunc)VectorObject_neg,        /*nb_negative*/
    0,                       /*nb_positive*/
    (unaryfunc)VectorObject_length,        /*nb_absolute*/
    (inquiry)VectorObject_nonzero,       /*nb_nonzero*/
    (unaryfunc)VectorObject_neg,           /*nb_invert*/
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

static PyTypeObject VectorObjectType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "VectorMath.Vector",       /*tp_name*/
    sizeof(VectorObject),      /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)VectorObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    VectorObject_cmp,          /*tp_compare*/
    0,                         /*tp_repr*/
    &VectorObject_as_number,    /*tp_as_number*/
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
    VectorObject_methods,      /* tp_methods */
    0,                         /* tp_members */
    VectorObject_getset,       /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)VectorObject_init, /* tp_init */
    0,                           /* tp_alloc */
    VectorObject_new            /* tp_new */
};

/*
 ***********************************************************************
 * 
 * Core functions of VectorMath library
 * 
 ***********************************************************************
 */
 
static PyObject *VectorObject_add(VectorObject *self, VectorObject *other) {
	
	VectorObject *result;    // Object being returned
	Vector *vr, *temp;  // C-vectors of a, b, and return value
	
	vr = vectorAdd(self->v, other->v);  // Add the two vectors and get the new one in return
	 
	if (!vr) {    // Ensure we recieved a vector 
		return NULL;
	}
	
	result    = (VectorObject *) VectorObject_new(&VectorObjectType, NULL); // Create new PyVector
	temp      = result->v;      // this whole process basically frees the vector currently allocated, and sets the pointer equal to the resultant vector
	result->v = vr;	            // It's faster than copying the data over because either way we have to free a vector. So this way we just toss a pointer around instead of 3 doubles
	
	vectorFree(temp);
	
	return (PyObject *) result;
}

static PyObject *VectorObject_sub(VectorObject *self, VectorObject *other) {
	
	VectorObject *result;    // Object being returned
	Vector *vr, *temp;  // C-vectors of a, b, and return value
	
	vr = vectorSub(self->v, other->v);  // Add the two vectors and get the new one in return
	 
	if (!vr) {    // Ensure we recieved a vector 
		return NULL;
	}
	
	result    = (VectorObject *) VectorObject_new(&VectorObjectType, NULL); // Create new PyVector
	temp      = result->v;
	result->v = vr;	
	
	vectorFree(temp);
	
	return (PyObject *) result;
}


static PyObject *VectorObject_length(VectorObject* self) {
	
    PyObject *result;
    double length;
   
    length = vectorLength(self->v);
    
    result = PyFloat_FromDouble(length);
    
    return result;
}

static PyObject *VectorObject_crossProduct(VectorObject *self, VectorObject *other) {
	
	VectorObject *retObject;
	Vector *vr, *temp;
	
	vr = crossProduct(self->v, other->v); // Calculate cross-product
	
	if (!vr) { // ensure we have a vector (High probability that we don't since cross product only works on 3D vectors)
		return NULL;
	}
	
	retObject = (VectorObject *) VectorObject_new(&VectorObjectType, NULL);
	temp      = retObject->v;
	retObject->v = vr;
	
	vectorFree(temp);
	
	return (PyObject *) retObject;
}

static PyObject *VectorObject_dotProduct(VectorObject *self, VectorObject *other) {
	
	PyObject *retObject;
	double result;
	
	result = dotProduct(self->v, other->v); // Calculate scalar-product
	
	retObject = PyFloat_FromDouble(result);
	
	return retObject;
}

static PyObject *VectorObject_angle(VectorObject *self, VectorObject *other) {
	
	PyObject *retObject;
	double result;
	
	result = angle(self->v, other->v); // Calculate scalar-product
	
	retObject = PyFloat_FromDouble(result);
	
	return retObject;
}


static PyObject *VectorObject_copy(VectorObject *self) {
	
	VectorObject *result = (VectorObject *) VectorObject_new(&VectorObjectType, NULL);
	Vector *temp, *v = vectorCopy(self->v);
	
	temp = result->v;
	result->v = v;
	
	free(temp);
	
	return (PyObject *) result;
}

static int VectorObject_cmp(VectorObject *self, VectorObject *other) {
	
	if (!PyObject_TypeCheck(other, &VectorObjectType)) { 
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

static int VectorObject_nonzero(VectorObject *self) {
	
	if (vectorLength(self->v)) { // If it has any sort of length, then it is not zero
		return 1;
	}
	
	return 0;  // Otherwise it is zero
}

static PyObject *VectorObject_mul(VectorObject *self, PyObject *other) {
	
	double k;
	VectorObject *ret;
	Vector *temp, *vr;
	
	if (PyObject_TypeCheck(other, &VectorObjectType)) {  // If the other object is a vector
	
		return VectorObject_dotProduct(self, other);     // make the crude assumption they meant dot product.
		 
	} else {
		
		k   = PyFloat_AsDouble(other);          // Otherwise, we assume it's either a float or an int
		ret = (VectorObject *) VectorObject_new(&VectorObjectType, NULL);  // Create new pyvector object
		
		vr     = vectorMul(self->v, k);       // multiply the vector by constant k
		temp   = ret->v;              // set temp equal to currently allocated vector
		ret->v = vr;                  // set the vector pointer to the resulting vector
		
		vectorFree(temp);             // free the previously allocated vector
		
	    return (PyObject *) ret;
	}
}

static PyObject *VectorObject_div(VectorObject *self, PyObject *other) {
	
	double k;
	VectorObject *ret;
	Vector *temp, *vr;
	
	if (PyObject_TypeCheck(other, &VectorObjectType)) {  // If the other object is a vector
	
		return VectorObject_dotProduct(self, other);     // make the crude assumption they meant dot product.
		 
	} else {
		
		k   = PyFloat_AsDouble(other);          // Otherwise, we assume it's either a float or an int
		ret = (VectorObject *) VectorObject_new(&VectorObjectType, NULL);  // Create new pyvector object
		
		vr     = vectorDiv(self->v, k);       // multiply the vector by constant k
		temp   = ret->v;              // set temp equal to currently allocated vector
		ret->v = vr;                  // set the vector pointer to the resulting vector
		
		vectorFree(temp);             // free the previously allocated vector
		
	    return (PyObject *) ret;
	}
}

static PyObject *VectorObject_neg(VectorObject *self) {
	
	VectorObject *ret;
	Vector *inverse, *temp;
	
	ret = VectorObject_new(&VectorObjectType, NULL); // create PyVector object
	
	inverse = vectorSub(ret->v, self->v);  // new vector is initialized to zero, so we subtract the vector to get it's inverse
	temp    = ret->v;
	ret->v  = inverse;
	
	vectorFree(temp);  // free allocated vector
	
	return (PyObject *) ret;
}

static PyObject *VectorObject_normalize(VectorObject *self) {
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
        
    if (PyType_Ready(&VectorObjectType) < 0) {     // Finalize type
		return;
	}
        
    Py_INCREF(&VectorObjectType);                  // Increment reference to pass to Python
    PyModule_AddObject(m, "Vector", (PyObject *)&VectorObjectType);  // Add to objects list so we can see it
}
