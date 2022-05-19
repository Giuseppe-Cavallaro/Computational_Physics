#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "phi4.h"

/* 
   allocates memory for phi, mom and hop arrays
   and calls hopping() 
*/
extern void initialize_fields();

/*
  reallocates memory of fields and calls hopping()
 */
extern void realloc_fields();

/* HOPPING_C */
extern void hopping(int (*h)[2*D]);

#endif
