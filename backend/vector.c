# include "vector.h"
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

DLLEX Vector *vectorNew(double *components, int demensions) {
	
	int componentSize = demensions*sizeof(double); // size of components in bytes
	Vector *v = malloc(sizeof(Vector));            // allocate memory
	
	v->demensions = demensions;                    // set number of demensions
    
    memcpy(v->components, components, componentSize);  // copy components to new vector
	
	return v;
}

DLLEX Vector *vectorCopy(Vector *a) {
	
	return vectorNew(a->components, a->demensions); // Just create a new vector, as it copies the components for us
	
}

DLLEX void vectorFree(Vector *v) {   // Just wraps free(). Had purpose when this supported infinite demension computations
	free(v);
}

DLLEX Vector *vectorAdd(Vector *a, Vector *b) {
	
	Vector *v;

	v = vectorCopy(a);    // create a vector copy of a, which will be returned 
	
	// add up matching components
	v->components[0] += b->components[0];
	v->components[1] += b->components[1];
	v->components[2] += b->components[2];
	
	return v;
}

DLLEX Vector *vectorSub(Vector *a, Vector *b) {
	
	Vector *v;
	
	v = vectorCopy(a);
	
	// subtract matching components
	v->components[0] -= b->components[0];
	v->components[1] -= b->components[1];
	v->components[2] -= b->components[2];
	
	return v;
}

DLLEX double vectorLength(Vector *a) {
	
	double sum;
	
	// Vector length is equal to (i^2 + j^2 + k^2)^(1/2)
	// add the squares of all the components
	sum = pow(a->components[0], 2.0);
	sum += pow(a->components[1], 2.0);
	sum += pow(a->components[2], 2.0);
	
	return sqrt(sum); // return the square root of the sum, which is the length
}

DLLEX int vectorCompare(Vector *a, Vector *b) {
	
	double l1, l2;
	
	l1 = vectorLength(a);  // Compute lengths of a and b
	l2 = vectorLength(b);
	
	if (l1 < l2) {        // Compare the length of each and return the appropriate code
		return -1;
	} else if (l1 == l2) {
		return 0;
	} else {
		return 1;
	}
}

DLLEX double dotProduct(Vector *a, Vector *b) { // scalar product
	
	int sum;
	
	// add up the products of the components for each vector
	sum = (a->components[0] * b->components[0]);
	sum += (a->components[1] * b->components[1]);
	sum += (a->components[2] * b->components[2]);
	
	return sum;
}

DLLEX Vector *crossProduct(Vector *a, Vector *b) { // vector product
	
	double components[3];
	
	// Use matrix method in order to compute the components of cross product
	components[2] = (a->components)[0]*(b->components)[1] - (a->components)[1]*(b->components)[0];
	components[0] = (a->components)[1]*(b->components)[2] - (a->components)[2]*(b->components)[1];
	components[1] = (a->components)[2]*(b->components)[0] - (a->components)[0]*(b->components)[2];
	
	return vectorNew(components, 3); // Return new vector using calculated components
}

DLLEX double angle(Vector *a, Vector *b) {
	
	double theta, dot, lenA, lenB;;
	
	// We know that a.b = |a||b|cos (theta), so by computing the 
	// dot-product, and the magnitudes of a and b, we can solve for cos (theta)
	// by rearranging to cos (theta) = a.b/(|a||b|)
	// and use acos to find theta
	
	dot  = dotProduct(a, b);
	lenA = vectorLength(a); 
	lenB = vectorLength(b);
	
	theta = acos(dot/(lenA*lenB));
	
	return theta * 180.0 / PI;
}


#ifdef BASIC_TIME_TEST // Little loop I made up for comparing pure C vs. Python vs. ctypes vs. extension

# include <time.h>

void main() {
	int i;
	int start;
	double c1[3] = {22.5, 12.8, 80.2};
	double c2[3] = {81.0, 5.0, 56.1};
	Vector *a = vectorNew(c1, 3);
	Vector *b = vectorNew(c2, 3);
	start = clock();
	for (i=0; i<100000; i++) {
		vectorAdd(a, b);
		vectorSub(a, b);
		vectorCompare(a, b);
		angle(a, b);
		vectorLength(a);
		vectorLength(b);
		dotProduct(a, b);
		crossProduct(a, b);
	}
	printf("%f ", (float)(clock()-start)/CLOCKS_PER_SEC);
	getchar();
}

#endif
