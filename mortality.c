/***********************************************************************************/
/* Program: mortality.c
   By: Brad Duthie                                             
   Description: Kills off individuals from both the ID and OFF matrices
   Compile: gcc mortality.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void mortality(double **ID, double **OFF, int Nloci, int Liv, int l, int M, 
    double **SCAPE, double Beta1, int DOM, int loadstart, int load, double cost,
    double hval){

    int i, k, a1, a2;
    double survival;

    /* First kill the older generation individuals in Liv */
    for(i=0; i<Liv; i++){ /* If individuals are too old (over M -- usually zero), die */
        if(ID[i][4]>=(M+1)){ /* They now have no location and -1 age */
            ID[i][2] = -1;
            ID[i][3] = -1;
            ID[i][4] = -1; 
        }
        if(ID[i][2] < 0){    /* And, for added measure, ensure they stay dead below */
            ID[i][4] = -1;
        }
    }    

    for(i=0; i<l; i++){ /* Next figure out if anyone dies in the new generation, OFF */
        survival = 1.0; /* Survival starts at 1, reduced by del recessives */
        for(k=0; k<load; k++){ /* For each inbreeding load loci */
            a1 = OFF[i][((4*k)+loadstart)];   /* Get the allele at position 1   */
            a2 = OFF[i][((4*k)+loadstart)+1]; /* Get the allele at position 2   */
            if(a1 == 0 && a2 == 0){           /* If both deleterious recessive  */
                survival *= (1 - Beta1);  /* Decrease survival probability  */
            }else{                            /* If both are not del recessive  */
                if(a1 == 0 || a2 == 0){   /* Is at least one del recessive? */
                    survival *= (1 - hval*Beta1);
                } /* If so, decrease survival probability (depends on hval) */
            }    
        } /* Value of survival is now the probability offspring survives (all loci) */
        OFF[i][2]  = survival; /* Put in column 2 for now (est. inbr. dep in sumstat)*/
        if(randunif()>survival){ /* Check to see if the individual dies */
            OFF[i][4] = -1;  /* If so, then add the -1 in column 4 */
        }
    } /* All juveniles should now have either survived, or died as juveniles */
}


