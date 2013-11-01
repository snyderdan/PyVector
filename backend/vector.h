# ifndef VECTORMATH_H
# define VECTORMATH_H

typedef double vectorComponent;

typedef struct _vector {
	int demensions;                // number of demensions; carried over from previous implementation allowing for greater than three demensional processing
	vectorComponent components[3]; // components of vector
} Vector;

Vector *vectorNew(vectorComponent *components, int demensions);   // create a vector based on the components and number of demensions specified
void   vectorFree(Vector *vector);         // free vector
Vector *vectorCopy(Vector *v);             // return copy of vector

# define PI 3.14159265359 

# define DLLEX __declspec(dllexport)    // declare for DLL compiling

DLLEX Vector __cdecl *vectorAdd(Vector *a, Vector *b);    // Adds vectors A and B
DLLEX Vector __cdecl *vectorSub(Vector *a, Vector *b);    // Subtracts vector b from vector a
DLLEX int    __cdecl vectorCompare(Vector *a, Vector *b); // Compares the magnitude of two vectors, returns: -1 if a < b; 0 if a == b; 1 if a > b
DLLEX double __cdecl vectorLength(Vector *v);             // Returns magnitude or length of the vector
DLLEX Vector __cdecl *vectorMul(Vector *a, double b);     // Multiply the vector by a scalar
DLLEX Vector __cdecl *vectorDiv(Vector *a, double b);     // Divide vector by scalar
DLLEX double __cdecl dotProduct(Vector *a, Vector *b);    // Calculates and returns the dot product of the two vector quantities; a and b remain unaltered
DLLEX Vector __cdecl *crossProduct(Vector *a, Vector *b); // Calculates the cross or scalar product of a * b and stores in a
DLLEX double __cdecl angle(Vector *a, Vector *b);   // Returns the angle between the two vectors

# endif
