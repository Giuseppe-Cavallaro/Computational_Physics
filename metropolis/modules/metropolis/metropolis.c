	
/*
 *  File metropolis.c
 *  Metropolis algorithm for the phi4 theory
 *
 *  double action(act_params_t *apars)
 *      This routine computes the action S[phi] for the global field phi in
 *      lattice.h and the parameters kappa and lambda from apars.
 *      S = Sum_x [ -2*kappa*sum_mu phi_x phi_{x+mu}+phi_x^2+lambda(phi_x^2-1)^2 ]
 */

#include <time.h>
#include "metropolis.h"

#define LBAR 50        /* length of loading bar*/
#define BAR_CHAR '='   /* character to print in loading bar */
#define FOLDER "../data/" /* folder for files containing observables data */
 
#define ANSI_COLOR_RED    "\x1b[31m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_RESET  "\x1b[0m"

char* attach_folder(char* filename, char* foldername) {
	char * path;
	if (filename == NULL) return NULL;
	path = (char *) malloc(1 + strlen(filename) + strlen(foldername) );
	strcpy(path, foldername);
	strcat(path, filename);
	return path;
}

double action(act_params_t *apars)
{
	int i, mu;
	double phin,S,phi2;
	double kappa =apars->kappa;
	double lambda=apars->lambda;

	S=0;

	/* loop over all sites */
	for (i=0;i<V;i++)
	{
		/*sum over neighbors in positive direction*/
		phin=0;
	    for (mu=0;mu<D;mu++) phin+=phi[hop[i][mu]];

		phi2=phi[i]*phi[i];
		S+=-2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
	}

	return S;
}

/* difference between actions calculated in phi prime (proposed for metropolis) and phi */
double deltaS(act_params_t *apars, metro_params_t *mpars, int i, double r) {
	double kappa  = apars->kappa;
	double lambda = apars->lambda;
	double delta  = mpars->delta;

	int mu;
	double muSum = 0;
	double p  = phi[i];
	double p2 = p*p;
	double p3 = p*p*p;
	double R  = delta * (r - 0.5);
	double R2 = R*R;
	double R3 = R*R*R;
	double R4 = R*R*R*R;

	for (mu = 0; mu < D; mu++) {
		muSum += phi[hop[i][mu]];
		muSum += phi[hop[i][mu+D]];
	}
	muSum *= -2 * R * kappa;

	return muSum + 2*p*R + R2 + lambda * (4*p3*R + 6*p2*R2 + 4*p*R3 - 4*p*R + R4 - 2*R2);
}

/* OBSERVABLES */

double M(){
	double sum = 0;
	int i;
	for (i = 0; i < V; i++) sum += phi[i];
	
	return sum;
}



/********************************************************************** 
 * metropolis()
 * Does ntherm+nsweep sweeps over the whole lattice of the Metropolis
 * algorithm. Measurement after each sweep>=ntherm, the averaged measured 
 * values are printed out in fixed intervals, controlled by naccu.
 **********************************************************************/

/* 
  * binX is the number of points to check with metropolis for every sweep:
  * it goes from 1 to V. If binX = V, metropolis checks every point on the lattice
  * for every sweep
*/
double metropolis(act_params_t *apars, metro_params_t *mpars, char *obs_file, int binX) {  
	double delta  = mpars->delta;
	double kappa  = apars->kappa;
	double lambda = apars->lambda;
	int Nsweep    = mpars->nsweep;
	int Ntherm    = mpars->ntherm;
	int binsize   = mpars->naccu; 
	
	int n, i;
	int nAcc  = 0; /* number of accepted values generated by metropolis */
	int nProp = 0; /* number of total proposed values by metropolis */
	double exp_dS, r;   

	int nM = 0;
	double tempM;
	double avgM = 0.0, avgMabs = 0.0, avgM2 = 0.0, avgM4 = 0.0; 

	FILE *f = fopen(attach_folder(obs_file, FOLDER), "w"); /* open file to write observable result */  

	time_t start, totStart;
	double seconds;
	int minutes;

	int j, cBar = 1;
	double ratio;
	unsigned char bar[LBAR + 1];  
	for (j = 0; j < LBAR; j++) bar[j] = ' ';
	bar[LBAR] = '\0';
	
	srand(time(NULL));
	r = random(); /* first generated number appears to be always the same */

	if (binX == V) {
		if (f == NULL) printf("Metropolis STARTED (binX = %d): " ANSI_COLOR_RED "\nNo file for printing observables" ANSI_COLOR_RESET "\n", binX);
		else printf("\nMetropolis STARTED (binX = %d): \nPrinting observables to %s\n", binX, obs_file);
	}
	else {
		if (f == NULL) printf("Metropolis STARTED (" ANSI_COLOR_YELLOW "binX = %d" ANSI_COLOR_RESET "): " ANSI_COLOR_RED "\nNo file for printing observables" ANSI_COLOR_RESET "\n", binX);
		else printf("\nMetropolis STARTED (" ANSI_COLOR_YELLOW "binX = %d" ANSI_COLOR_RESET ")\nPrinting observables to %s\n", binX, obs_file);
	}
	
	totStart = time(NULL);
	start = time(NULL);

	/* METROPOLIS THERMALIZATION (No observables calculation) */
	for (n = 0; n < Ntherm; n++) {
	    /* cycle through binX number of lattice points */
	    for (i = binX*n % V; i < binX*n % V + binX; i++){
			nProp++;
			r = random(); /* generate random number in (0,1) */
			exp_dS = exp( -deltaS(apars, mpars, i % V, r) ); /* calculate the exponential of -deltaS */

			if ((exp_dS >= 1) || (random() < exp_dS)) {
				phi[i % V] += delta * (r - 0.5); /* accept the proposed value for phi[i] */
				nAcc++;
			}
		}      

		/* calculating ETA */      		
		seconds = difftime(time(NULL), start);
		seconds = (seconds / (n+1)) * (Ntherm - (n+1));
		minutes = seconds / 60;
		seconds -= minutes * 60;
		
		/* loading bar */
		ratio = (double)(n+1) / (double)Ntherm;		
		if (ratio * LBAR >= cBar) {
			bar[cBar - 1] = BAR_CHAR;
			cBar++;
		}
		
		printf("\rTherm |%s| ETA: %02d:%02d %6.2f%%", bar, minutes, (int)seconds, ratio * 100); 
	}   
	printf("\n");

	/* reset bar */
	cBar = 1;
	for (j = 0; j < LBAR; j++) bar[j] = ' ';

	/* print file header */
	if (f != NULL) {
		fprintf(f, "%d %d %d %d %d %f %f %f\n", V, binsize, binX, Ntherm, Nsweep, delta, kappa, lambda); 
		fprintf(f, "n M absM M2 M4\n"); /* print name of observables */
	}

	start = time(NULL);

	/* METROPOLIS SWEEPING (Observables calculation) */
	for (n = 0; n < Nsweep; n++) {
		/* cycle through binX number of lattice points */
		for (i = binX*n % V; i < binX*n % V + binX; i++){
			nProp++;
			r = random(); /* generates random real number in (0,1) */
			exp_dS = exp( -deltaS(apars, mpars, i % V, r) ); /* calculates the exponential of -deltaS */
			
			if ((exp_dS >= 1) || (random() < exp_dS)) {
				phi[i % V] += delta * (r - 0.5); /* accept the proposed value for phi[i] */
				nAcc++;
			}
		}
		
		/* observables calculation */
		tempM    = M();
		avgM    += tempM;			
		avgMabs += fabs(tempM);		
		avgM2   += tempM * tempM;	
		avgM4   += tempM * tempM * tempM * tempM;
		
		if ((n+1) % binsize == 0) { /* true every binsize number of sweeps */
			nM++;
			avgM    /= binsize; 
			avgMabs /= binsize; 
			avgM2   /= binsize;
			avgM4   /= binsize;

			/* print observables to file */
			if (f != NULL) fprintf(f, "%d %e %e %e %e\n", nM, avgM, avgMabs, avgM2, avgM4);

			avgM    = 0.0;
			avgMabs = 0.0;
			avgM2   = 0.0;
			avgM4   = 0.0;
		}
		
		/* calculating ETA */
		seconds = difftime(time(NULL), start);
		seconds = (seconds / (n+1)) * (Nsweep - (n+1));
		minutes = seconds / 60;
		seconds -= minutes * 60;

		/* loading bar */
		ratio = (double)(n+1) / (double)Nsweep;		
		if (ratio * LBAR >= cBar) {
			bar[cBar - 1] = BAR_CHAR;
			cBar++;
		}
		
		printf("\rSweep |%s| ETA: %02d:%02d %6.2f%%", bar, minutes, (int)seconds, ratio * 100); 
	} 

	if (f != NULL) {
		fprintf(f, "pa %f\n", (double) nAcc / (double) nProp);
		fclose(f);
	}
	
	seconds = difftime(time(NULL), totStart);
	minutes = seconds / 60;
	seconds -= minutes * 60;

	printf("\nMetropolis COMPLETED in %02d:%02d\n\n", minutes, (int)seconds);

	return (double) nAcc / (double) nProp;
}
