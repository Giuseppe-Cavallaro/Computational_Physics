#include "math.h"
#include "my_rand.h"
#include "ranlxd.h"

#define PI 3.14159265

/* random number between 0 and 1 */ 
double random(){
	return ((double) rand() / (RAND_MAX));
}

/* fill a vector with gaussian distributed variables */
void gauss_rand(double *v, int dim) {
	int i;
	double temp[2];
	
	for (i = 0; i < dim; i += 2) {
		ranlxd(temp, 2);
		v[i]   = sqrt(-2*log(1-temp[0])) * sin(2*PI*(1-temp[1]));
		v[i+1] = sqrt(-2*log(1-temp[0])) * cos(2*PI*(1-temp[1]));
	}
	if (dim % 2 == 1) {
		ranlxd(temp, 2);
		v[dim-1] = sqrt(-2*log(1-temp[0])) * sin(2*PI*(1-temp[1]));
	}
}
