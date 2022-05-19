/*
 *   File check1.c
 *
 *   Checks reversibility of the algorithm
 */

#define CONTROL
#include "phi4.h"
#include "metropolis.h"
#include "hmc.h"

/*
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in phi4.h
 *  seed for initializing the random numbers
 */

static int seed = 147;

act_params_t apars;
hmc_params_t pars;

int main(int argc, char* argv[])
{
	int i;
	double *phi_i;
	double *mom_i;
	double diffphi = 0, diffmom = 0;
	double H_i, H_f;

	/* initialize random number generator */
	rlxd_init(1,seed);

	L = LATTICE_L;
	initialize_fields();
	realloc_fields();

	phi_i = malloc(sizeof(double) * V);
	mom_i = malloc(sizeof(double) * V);

	apars.kappa  = 0.15;
	apars.lambda = 1.145;

	pars.nstep   = 10;
	pars.tlength = 1;

	printf("PARAMETERS\n");
	printf("L           %i\n", L);
	printf("D           %i\n", D);
	printf("kappa       %f\n", apars.kappa);
	printf("lambda      %f\n", apars.lambda);
	printf("tlength     %f\n", pars.tlength);
	printf("nstep       %i\n", pars.nstep);
	printf("END PARAMETERS\n");
	printf("\n");

	ranlxd(phi, V);
	gauss_rand(mom, V);

	for (i = 0; i < V; i++) {
		phi_i[i] = phi[i];
		mom_i[i] = mom[i];
	}

	H_i = hamilton(&apars);
	printf("Initial H = %f\n", H_i);

	leapFrog(&apars, &pars);

	printf("Middle H = %f\n", hamilton(&apars));

	for (i = 0; i < V; i++) mom[i] *= -1;

	leapFrog(&apars, &pars);

	H_f = hamilton(&apars);
	printf("Final H = %f\n", H_f);

	for (i = 0; i < V; i++) {
		diffphi += fabs(phi[i] - phi_i[i]);
		diffmom += fabs(mom[i] - mom_i[i]);
	}

	diffphi /= V;
	diffmom /= V;

	printf("\nPhi diff = %e\n", diffphi);
	printf("Mom diff = %e\n", diffmom);
	printf("H diff = %e\n", H_f - H_i);

	return 0;
}
