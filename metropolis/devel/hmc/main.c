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
static int minL, maxL, dL;
static double min_kappa, max_kappa, d_kappa;

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

   get_val(fp, "lambda",      "%lf", &act_params.lambda);
 
   get_val(fp, "min_kappa",   "%lf", &min_kappa);
   get_val(fp, "max_kappa",   "%lf", &max_kappa);
   get_val(fp, "dkappa",      "%lf", &d_kappa);
   
   get_val(fp, "min_L",       "%i" , &minL);
   get_val(fp, "max_L",       "%i" , &maxL);
   get_val(fp, "dL",          "%i" , &dL);

   get_val(fp, "tlength",     "%lf", &hmc_params.tlength);
   get_val(fp, "ntherm",      "%i" , &hmc_params.ntherm);
   get_val(fp, "ntraj",       "%i" , &hmc_params.ntraj);
   get_val(fp, "naccu",       "%i" , &hmc_params.naccu);
   get_val(fp, "seed",        "%i" , &seed);
   			
   printf("PARAMETERS\n");
   printf("DIM             %i\n",     D                    );
   printf("L in range      %i, %i\n", minL, maxL           );
   printf("L increment     %i\n",     dL                   );
   printf("kappa in range  %f, %f\n", min_kappa, max_kappa );
   printf("kappa increment %f\n",     d_kappa              );
   printf("lambda          %f\n",     act_params.lambda    );
   printf("tlength         %f\n",     hmc_params.tlength   );
   printf("ntherm          %i\n",     hmc_params.ntherm    );
   printf("ntraj           %i\n",     hmc_params.ntraj     );
   printf("naccu           %i\n",     hmc_params.naccu     );
   printf("END PARAMETERS\n");

   return 0;
}

/* get the number of step of LeapFrog() to have acc rate ~ 0.9 */
static int get_nstep(double c) {
	return (int) ceil( c * sqrt( sqrt( pow(L, D) ) ) );
}

int main(int argc, char* argv[]) {
	double c = 3.0;
	double acc;
	double *pa = NULL;
	double k;
	int i = 0, N_pa = 0;
	FILE *f;

	time_t start;
	int minutes;
	double seconds;
	
	if (argc < 3) 
	{
		fprintf(stderr,"Number of arguments not correct\n");
		fprintf(stderr,"Usage: %s <infile> <output file>\n",argv[0]);
		exit(1);
	}

	/* set field pointers to NULL */
	initialize_fields();
	
	/* get the parameters from the input file */
    read_input(argv[1]);
	printf("\n");
	
    /* initialize random number generator */
    rlxd_init(1, seed);
	srand(time(NULL));

	/* calculates N_pa */
	for (L = minL; L <= maxL; L += dL) {
		for (k = min_kappa; k <= max_kappa; k += d_kappa) {
			N_pa++;
		}
	}
	pa = malloc(sizeof(double) * N_pa);	
	
	/* open output file */
	f = fopen(attach_folder(argv[2], FOLDER), "w");	
	if (f == NULL) {
		printf("error opening file\n");
		return 1;
	}

	printf("Input file: %s\n", argv[1]);
	printf("Output file: %s%s\n", FOLDER, argv[2]);
	printf("\n");

	/* print info about program on the first two lines  */
	fprintf(f, "D,lambda,tlength,ntherm,ntraj,naccu,c\n");
	fprintf(f, "%i,%e,%e,%i,%i,%i,%e\n",
			D,
			act_params.lambda,
			hmc_params.tlength,
			hmc_params.ntherm,
			hmc_params.ntraj,
			hmc_params.naccu,
			c);
	fprintf(f, "\n"); /* leave a blank line in file */
	fprintf(f, "L,k,m,mAbs,m2,m2_V,m4\n"); /* print observables header */

	start = time(NULL);
	
	/* does HMC for every L and k in the range given */
	for (L = minL; L <= maxL; L += dL) {
		printf("HMC for L = %i...\n", L);
		hmc_params.nstep = get_nstep(c);
		realloc_fields();
		for (k = min_kappa; k <= max_kappa; k += d_kappa) {   
			act_params.kappa = k;
			ranlxd(phi, V);
			acc = hmc(&act_params, &hmc_params, f);
			pa[i] = acc;
			i++;
			printf("   k = %f - acc = %f\n", k, acc);
		}
		printf("\n");
	}

	/* print probability of acceptance to file */
	fprintf(f, "\n");
    fprintf(f, "L,k,pa\n");
	i = 0;
	for (L = minL; L <= maxL; L += dL) {
		for (k = min_kappa; k <= max_kappa; k += d_kappa) {
			fprintf(f, "%i,%e,%e\n", L, k, pa[i]);
			i++;
		}
	} 

	seconds = difftime(time(NULL), start);
	minutes = seconds / 60;
	seconds -= minutes * 60;

	printf("HMC completed in %02d:%02d\n\n", minutes, (int)seconds);
	
	fclose(f);
	free(phi);
	free(mom);
	free(hop);
	free(pa);

	return 0;
}

