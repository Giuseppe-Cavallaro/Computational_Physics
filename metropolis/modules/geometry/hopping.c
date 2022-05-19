/*********************************************************
 * 
 *  File hopping.c
 *
 *  Initialization of the hopping field for a D dimensional
 *  lattice of size V = L**D
 *  The index of a point with coordinates (n_0,n_1,..,n_{D-1})
 *  is i = sum_k n_k L**k
 *  The index of its neighbor in positive direction nu 
 *  is hop[i][mu]
 *  In negative direction it is hop[i][D+mu]
 *
 ********************************************************/

#include "phi4.h"

void initialize_fields() {
	hop     = NULL;
	phi     = NULL;
	phi_old = NULL;
	mom     = NULL;	
}

void realloc_fields() {
	V = pow(L, D);
	hop = realloc(hop, sizeof(int[2*D]) * V); 
	phi = realloc(phi, sizeof(double) * V);
	phi_old = realloc(phi_old, sizeof(double) * V);
	mom = realloc(mom, sizeof(double) * V);	
	hopping(hop);	
}

void hopping(int (*hop)[2*D])
{
   int x, y, Lk;
   int xk, k, dxk;

   /* go through all the points*/
   for (x = 0; x < V; x++)
   {
      Lk = V;
      y = x;

      /* go through the components k*/
      for (k = D-1; k >= 0; k--)
      {
         Lk /= L;                        /* pow(l,k)      */ 
         xk = y / Lk;                    /* kth component */
         y  = y - xk * Lk;               /* y<-y%Lk       */

         /* forward */
         if (xk < L-1) 
            dxk = Lk;
         else
            dxk = Lk * (1-L);
         hop[x][k] = x + dxk;

         /* backward */
         if (xk > 0)   
            dxk = -Lk;
         else
            dxk= Lk * (L-1);
         hop[x][k+D] = x + dxk;
      }
   }
} /* hopping */

