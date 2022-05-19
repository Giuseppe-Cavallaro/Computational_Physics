
/* 
 *   File check4.c
 *
 *   <exp(-deltaH)> = 1 in theory  
*/ 

#define CONTROL
#include "phi4.h"
#include "metropolis.h"
#include "hmc.h"
#include <time.h>

#define FOLDER "../data/"

/*  
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in phi4.h
 *  seed for initializing the random numbers
 */

static hmc_params_t hmc_params;
static act_params_t act_params;
static int seed;

static int get_val(FILE* fp, char *str, char* fmt,  void* val)
{
   char c[128];

   if(1!=fscanf(fp,"%s",c))
   {
      fprintf(stderr,"Error reading input file at %s\n",str);
      exit(1);
   }

   if(strcmp(str,c)!=0)
   {
      fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
      exit(1);
   }

   if(1!=fscanf(fp,fmt,val))
   {
      fprintf(stderr,"Error reading input file at %s\n",str);
      fprintf(stderr,"Cannot read value format %s\n",fmt);
      exit(1);
   }

   return 0;
}


static int read_input(char *input)
{
   FILE* fp;

   fp=fopen(input,"r");

   if (fp==NULL) 
   {
      fprintf(stderr, "Cannot open input file %s \n",input);
      exit(1);
   }

   get_val(fp, "kappa",       "%lf", &act_params.kappa);
   get_val(fp, "lambda",      "%lf", &act_params.lambda);
   get_val(fp, "L",           "%i" , &L);
   get_val(fp, "nstep",       "%i" , &hmc_params.nstep);
   get_val(fp, "tlength",     "%lf", &hmc_params.tlength);
   get_val(fp, "ntherm",      "%i" , &hmc_params.ntherm);
   get_val(fp, "ntraj",       "%i" , &hmc_params.ntraj);
   get_val(fp, "naccu",       "%i" , &hmc_params.naccu);
   get_val(fp, "seed",        "%i" , &seed);

   printf("PARAMETERS\n");
   printf("L              %i\n", L                   );
   printf("DIM            %i\n", D                   );
   printf("kappa          %f\n", act_params.kappa    );
   printf("lambda         %f\n", act_params.lambda   );
   printf("nstep          %i\n", hmc_params.nstep    );
   printf("tlength        %f\n", hmc_params.tlength  );
   printf("ntherm         %i\n", hmc_params.ntherm   );
   printf("ntraj          %i\n", hmc_params.ntraj    );
   printf("naccu          %i\n", hmc_params.naccu    );
   printf("END PARAMETERS\n");

   return 0;
}

double hmc_check4(act_params_t *apars, hmc_params_t *pars, FILE *f) {
	double ntherm = pars->ntherm;
	double ntraj  = pars->ntraj;
	int    naccu  = pars->naccu;
	
	int nProp = 0, nAcc = 0;
	int n, i;
	double *phi_old = malloc(sizeof(double) * V);
	double H_i, deltaH;

	double exp_deltaH_avg = 0.0;
	
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
		exp_deltaH_avg += exp(-deltaH);
		
		if ((deltaH >= 0) && (random() > exp(-deltaH))) {
			for (i = 0; i < V; i++) phi[i] = phi_old[i];
		}
		else nAcc++;

		if ((n+1) % naccu == 0) {
			exp_deltaH_avg /= naccu;
			if (f != NULL) fprintf(f, "%e\n", exp_deltaH_avg);
			exp_deltaH_avg = 0.0;
		}
	}
	
	return (double)nAcc / (double)nProp;
}

static void print_file_info(FILE *f, act_params_t apars, hmc_params_t pars) {
	fprintf(f, "D,L,kappa,lambda,nstep,tlength,ntherm,ntraj,naccu\n");
	fprintf(f, "%i,%i,%e,%e,%i,%e,%i,%i,%i\n",
			D, L, apars.kappa, apars.lambda, pars.nstep,
			pars.tlength, pars.ntherm, pars.ntraj, pars.naccu);
	fprintf(f, "\n");
}

int main(int argc, char* argv[]) {
	double acc;
	char *filename = "check4.csv";
	FILE *f;
	
	if ((argc < 2) || (argc > 3)) 
	{
		fprintf(stderr,"Number of arguments not correct\n");
		fprintf(stderr,"Usage: %s <infile> <output file>(optional)\n",argv[0]);
		printf("Input file format must be:\n");
		printf("kappa    \n");
		printf("lambda   \n");
		printf("L        \n");
		printf("nstep    \n");
		printf("tlength  \n");
		printf("ntherm   \n");
		printf("ntraj    \n");
		printf("naccu    \n");
		printf("seed     \n");
		exit(1);
	}
	else if (argc == 3) filename = argv[2]; 
	
	/* get the parameters from the input file */
    read_input(argv[1]);
	printf("\n");
	
    /* initialize random number generator */
    rlxd_init(1,seed);
	srand(time(NULL));

    initialize_fields();
	ranlxd(phi, V);
	
	f = fopen(attach_folder(filename, FOLDER), "w");	
	if (f == NULL) {
		printf("error opening file\n");
		return 1;
	}

	print_file_info(f, act_params, hmc_params);
	
	fprintf(f, "exp_deltaH\n");
	acc = hmc_check4(&act_params, &hmc_params, f);
    fclose(f);
	printf("acc = %f\n", acc);
    printf("\n");
	return 0;
}

