#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "phi4.h"

/* METROPOLIS_C */
extern char* attach_folder(char* filename, char* foldername);
extern double action(act_params_t *apars);
extern double M();
extern double metropolis(act_params_t *apars, metro_params_t *mpars, char *obs_file, int binX);

#endif
