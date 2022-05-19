#include "hmc.h"

/* F(x) = dS(phi)/dphi(x) */
double F(int x, act_params_t *apars) {
	double kappa  = apars->kappa;
	double lambda = apars->lambda;
	int mu;
	double muSum = 0.0;

	double p  = phi[x];
	double p2 = p*p;
	
	for (mu = 0; mu < D; mu++) {
		muSum += phi[hop[x][mu]];
		muSum += phi[hop[x][mu+D]];
	}

	return -2*kappa*muSum + 2*p + 4*lambda*p*(p2 - 1.0);
	
}

/* return the Hamiltonian with the action as the potential */
double hamilton(act_params_t *apars) {
	int i;
	double sumSqr = 0.0;
	
	for (i = 0; i < V; i++){
		sumSqr += mom[i]*mom[i]; 
	}
	
	return 0.5*sumSqr + action(apars);
}

void move_mom(double eps, act_params_t *apars) {
    int x;

	for (x = 0; x < V; x++) {
		mom[x] -= eps * F(x, apars);
	}
}

void move_phi(double eps) {
    int x;

	for (x = 0; x < V; x++) {
		phi[x] += eps * mom[x];
	}
}

void leapFrog(act_params_t *apars, hmc_params_t *pars) {
	int Nstep = pars->nstep;
	double dt = pars->tlength / Nstep;

	int i;

	/* apply Leap Frog Nstep times */	
	move_mom(dt / 2, apars);
	move_phi(dt);
	for (i = 0; i < Nstep - 1; i++) {
		move_mom(dt, apars);
		move_phi(dt);
	}
    move_mom(dt / 2, apars);
}

double hmc(act_params_t *apars, hmc_params_t *pars, FILE *f) {
	double ntherm = pars->ntherm;
	double ntraj  = pars->ntraj;    
	int    naccu  = pars->naccu;
	double k      = apars->kappa;
	
	int nProp = 0, nAcc = 0;
	int n, i;
	double H_i, deltaH;

	/* variable to store sums of observables */
	double tempM;
	double m    = 0.0;
	double mAbs = 0.0;
	double m2   = 0.0;
	double m2_V = 0.0;
	double m4   = 0.0;
	
	/* THERMALIZATION */
	for (n = 0; n < ntherm; n++) {
		nProp++;
		gauss_rand(mom, V);
		for (i = 0; i < V; i++) phi_old[i] = phi[i];
		H_i = hamilton(apars);
		leapFrog(apars, pars);
		deltaH = hamilton(apars) - H_i;
		if ((deltaH >= 0) && (random() > exp(-deltaH))) {
			for (i = 0; i < V; i++) phi[i] = phi_old[i];
		}
		else nAcc++;
	}

	/* TRAJECTORIES */
	for (n = 0; n < ntraj; n++) {
		nProp++;
		gauss_rand(mom, V);
		for (i = 0; i < V; i++) phi_old[i] = phi[i];
		H_i = hamilton(apars);
		leapFrog(apars, pars);
		deltaH = hamilton(apars) - H_i;
		if ((deltaH >= 0) && (random() > exp(-deltaH))) {
			for (i = 0; i < V; i++) phi[i] = phi_old[i];
		}
		else nAcc++;

		/* calculate observables */
		tempM = M();
		m    += tempM;
		mAbs += fabs(tempM);
		m2   += tempM*tempM;
		m2_V += tempM*tempM / V;
		m4   += tempM*tempM*tempM*tempM;
		
		if ((n+1) % naccu == 0) {
			m    /= naccu;
			mAbs /= naccu;
			m2   /= naccu;
			m2_V /= naccu;
			m4   /= naccu;
			if (f != NULL)
				fprintf(f, "%i,%e,%e,%e,%e,%e,%e\n", L, k, m, mAbs, m2, m2_V, m4);
			m    = 0.0;
			mAbs = 0.0;
			m2   = 0.0;
			m2_V = 0.0;
			m4   = 0.0;
		}
	}

	return (double)nAcc / (double)nProp;
}
