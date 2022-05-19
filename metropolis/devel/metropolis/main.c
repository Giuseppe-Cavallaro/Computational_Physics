
/* 
 *   File check1.c
 *
 *   Contains the main program and a few other routines from which the
 *   code to simulate the phi**4 theory can be built. Routines for reading 
 *   in the main parameters of the action and the algorithm are provided.
 *
 *   static int get_val(FILE* fp, char *str, char* fmt,  void* val)
 *      Routine which reads one line from the input file.
 *      Format of the lines is <keyword> <value>.
 *      Checks if the keyword in string str matches,
 *      then gets the value according to the format in fmt
 *
 *   static int read_input(char *input)
 *      Parses the input file (format as specified in get_val)
 *      and prints the parameters onto the screen. Currently
 *      it reads the basic values for the action and also for the 
 *      future Metropolis and the seed of the random number generator.
 */ 

#define CONTROL
#include "phi4.h"
#include "metropolis.h"

/*  
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in phi4.h
 *  seed for initializing the random numbers
 */

static metro_params_t metro_params;
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

   get_val(fp, "kappa",       "%lf",&act_params.kappa    );
   get_val(fp, "lambda",      "%lf",&act_params.lambda   );
   get_val(fp, "ntherm",      "%i" ,&metro_params.ntherm );
   get_val(fp, "nsweep",      "%i" ,&metro_params.nsweep );
   get_val(fp, "delta",       "%lf",&metro_params.delta  );
   get_val(fp, "seed",        "%i" ,&seed                );
   get_val(fp, "naccu",       "%i" ,&metro_params.naccu  );

   printf("PARAMETERS\n");
   printf("L              %i\n", L                   );
   printf("DIM            %i\n", D                   );
   printf("kappa          %f\n", act_params.kappa    );
   printf("lambda         %f\n", act_params.lambda   );
   printf("ntherm         %i\n", metro_params.ntherm );
   printf("nsweep         %i\n", metro_params.nsweep );
   printf("delta          %f\n", metro_params.delta  );
   printf("naccu          %i\n", metro_params.naccu  );
   printf("END PARAMETERS\n");

   return 0;
}

int main(int argc, char* argv[])
{
   double act,acc;
   char *obs_file = NULL; 
   int binX = LATTICE_V;
   L = LATTICE_L;
   V = LATTICE_V;
   
   if (argc < 2) 
   {
	  fprintf(stderr,"Number of arguments not correct\n");
	  fprintf(stderr,"Usage: %s <infile> \n",argv[0]);
	  exit(1);
   }
   switch (argc) {
      case 3:
	      obs_file = argv[2]; break;
      case 4: {
         obs_file = argv[2];
         binX = atoi(argv[3]);
         break;
      }
   }

   /* get the parameters from the input file */
   read_input(argv[1]);

   /* initialize random number generator */
   rlxd_init(1,seed);

   /* initialize the nearest neighbor field */
   hop = malloc(sizeof(int[2*D]) * V);
   hopping(hop);

   /* initialize phi field */
   phi = malloc(sizeof(double) * V);
   ranlxd(phi,V);
   
   /* compute the action */
   act = action(&act_params);
   printf("ACTION = %17.7e\n\n", act);

   /* do the updating */
   acc = metropolis(&act_params, &metro_params, obs_file, binX);
   printf("ACCRATE %f\n",acc);

   /* compute the action */
   act = action(&act_params);
   printf("ACTION = %17.7e\n", act);

   return 0;
}

