#ifndef LATTICE_H
#define LATTICE_H

/* dimension of the lattice */
#define D 3
/* spatial extend of the lattice: default value */
#define LATTICE_L 4
/* lattice volume, needs to be adjusted according to number of dimensions: default value */
#define LATTICE_V (LATTICE_L*LATTICE_L*LATTICE_L)

#ifdef CONTROL 
#define EXTERN 
#undef CONTROL
#else
#define EXTERN extern
#endif

EXTERN int L;
EXTERN int V;

EXTERN double *phi;
EXTERN double *phi_old;
EXTERN int    (*hop)[2*D];
EXTERN double *mom;

#endif
