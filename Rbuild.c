/***********************************************************************************/
/* Program: Rbuild.c
   By: Brad Duthie                                             
   Description: Uses kinship matrix from previous generation(s) to build new kinships.
   This is called Rbuild.c, but it really is calculating kinships, not relatedness

   Compile: gcc Rbuild.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"

void Rbuild(double **ID, double **Rmat, double **Rmof, int Liv, int OLiv, int M){

    int i, j; 
    int MomP, DadP, Om, Od;
    double rMmuP, rDmuP, rval, w, x, y, z;

    /* =====================================================================*/
    /* Below finds the correct kinship value between all individuals except */
    /* immigrants and their children                                        */
    /* =====================================================================*/
    for(i=0; i<Liv; i++){
        /* Find Mom position -- -1 if immigrant */
        if(ID[i][5] < 0){
            MomP = -1;
        }else{
            MomP = 0; /* printf("\t\t\t%f\n",Rmat[0][0]); */
            while(ID[i][5]!=Rmat[MomP][0]){
                MomP++;
            }
        }
        /* Find Dad position -- -1 if immigrant */
        if(ID[i][6] < 0){
            DadP = -1; 
        }else{
            DadP = 0; /* Below finds the matrix position of Dad */
            while(ID[i][6]!=Rmat[DadP][0]){
                DadP++;
            }
        } 

        if(MomP < 0 || DadP < 0){
            Rmof[i][i+1] = 1; 
        }else{
            Rmof[i][i+1] = 1 + Rmat[MomP][DadP+1]; 
        }

        for(j=0; j<i; j++){ /* Looking for R for each offspring... */
            /* Find r between i Mom and j Mom (x) and Dad (y) */
            /* Establish who the parents are =======================*/
            if(ID[j][5] < 0){
                Om = -1;
            }else{
                Om = 0; /* Other mom's position */
                while(ID[j][5]!=Rmat[Om][0]){
                    Om++;
                } /* Keeps looking up until hits Mom's Rmat pos */
            }
            if(ID[j][6] < 0){
                Od = -1;
            }else{
                Od = 0; /* Other dad's position */
                while(ID[j][6]!=Rmat[Od][0]){
                    Od++; 
                } /* Keeps looking up until hits Dad's Rmat pos */
            }
            if(ID[i][4] >= 0 && ID[i][4] <= M && ID[j][4] >= 0 && ID[j][4] <= M){
                if(Om < 0 || MomP < 0){ /* If one Mom Imm */
                    w = 0; 
                }else{ /* Must be unrelated (not on the table) */
                    if(MomP == Om){
                        w = 0.5 * Rmat[MomP][Om+1];
                    }else{
                        w = Rmat[MomP][Om+1];
                    }
                }
                if(Od < 0 || MomP < 0){ /* Same with Mom1 and Dad2 */
                    y = 0;
                }else{
                    y = Rmat[MomP][Od+1]; 
                }
                if(Om < 0 || DadP < 0){
                    x = 0;
                }else{
                    x = Rmat[DadP][Om+1]; 
                }
                if(Od < 0 || DadP < 0){
                    z = 0;
                }else{
                    if(DadP == Od){
                        z = 0.5 * Rmat[DadP][Od+1];
                    }else{
                        z = Rmat[DadP][Od+1]; 
                    }
                }
                
                /* Now do the actual calculation */
                rMmuP = 0.5 * (w + y);
                rDmuP = 0.5 * (x + z);
                rval  = 0.5 * (rMmuP + rDmuP); 
            }else{ /* Other rvals between individuals will not be further used */
                rval  = 0.0; /* This includes at least one dead individual, now */
            } /* unable to again reproduce, so rval not needed for future offspring */
            Rmof[i][j+1] = rval;
            Rmof[j][i+1] = rval;
        }
    }

    /* Immigrants are only related to their descendents (par-off kinship = 0.25) */
    /* Below sets the parent-offspring kinship to 0.25 for children of immigrants */
    for(i=0; i<Liv; i++){
        for(j=0; j<Liv; j++){
            if(ID[i][5] == ID[j][0]){
                if(ID[j][5]==-1){
                    Rmof[i][j+1] = 0.25;
                    Rmof[j][i+1] = 0.25;
                } 
            }
            if(ID[i][6] == ID[j][0]){
                if(ID[j][5]==-1){
                    Rmof[i][j+1] = 0.25; 
                    Rmof[j][i+1] = 0.25; 
                }
            }
        }
    } 
}


