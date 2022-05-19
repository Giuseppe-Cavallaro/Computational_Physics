#ifndef PHI4_H
#define PHI4_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ranlxd.h"
#include "my_rand.h"
#include "lattice.h"
#include "geometry.h"

/*data structure to store all the parameters of the algorithm*/
typedef struct {
	double delta;  
	int ntherm;     /*number of thermalization sweeps*/
	int nsweep;     /*number of sweeps after thermalization*/
	int naccu;      /*binsize for printing out the measurements*/
} metro_params_t;

/*data structure to store all the parameters of the action*/
typedef struct {
	double kappa;
	double lambda;
} act_params_t;

/* data structure to store all the parameters of the hybrid montecarlo */
typedef struct {
	int nstep;     /* N zero */
	double tlength /* tau zero */;
	int ntherm;    /* # of thermalization step */
	int ntraj;     /* # of trajectories after thermalization */
	int naccu;     /* used for binning */
} hmc_params_t;

#endif
