#ifndef HYBRID_H
#define HYBRID_H

#include "phi4.h"
#include "metropolis.h"

extern double hamilton(act_params_t *apars);

/* mom -> mom - eps * dS/dphi(phi) */ 
extern void move_mom(double eps, act_params_t *apars);

/* phi -> phi + eps * mom */ 
extern void move_phi(double eps);

/* does one leapfrog trajectory */
extern void leapFrog(act_params_t *apars, hmc_params_t *pars);

/* implement hmc algorithm */
extern double hmc(act_params_t *apars, hmc_params_t *pars, FILE *f);

#endif
