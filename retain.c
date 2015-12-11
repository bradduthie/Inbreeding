/***********************************************************************************/
/* Program: retain.c
   By: Brad Duthie                                             
   Description: Marks individuals that are parents (so they aren't removed), 
        needed for Rmat calcs.
   Compile: gcc retain.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void retain(double **ID, double **OFF, int Liv, int l, int M){

    int i, j, k, h, g;
    /* Check if dead/too old to breed individuals are parents */    
    /*If they are, need to be retained for use in Rmat */
    h = 0; /* Checks to see if someone no longer parent of living individuals */
    do{  /* If they aren't, we can get rid of them permanently (-1) */
        g = 0;        
        for(i=0; i<Liv; i++){ /* Go through the entire ID array (adults) */
            if(ID[i][4] == -1){ /* If an individual has died, see if parent */
                k = 0; /* If have living child in ID, add to k */
                for(j=0; j<Liv; j++){
                    if(ID[i][0]==ID[j][5] && ID[j][4]>=0){
                        k++;
                    }
                    if(ID[i][0]==ID[j][6] && ID[j][4]>=0){
                        k++;
                    }
                } /* Or, if have living child in ID, add to k */
                for(j=0; j<l; j++){
                    if(ID[i][0]==OFF[j][5] && OFF[j][4]>=0){
                        k++;
                    }
                    if(ID[i][0]==OFF[j][6] && OFF[j][4]>=0){
                        k++;    
                    }
                }
                if(k == 0){ /* If they have no living children */
                    ID[i][4] = -2; /* Make them really, really dead */
                    ID[i][5] = -1; /* Their Dad is irrelevant too */
                    ID[i][6] = -1; /* And so is their Mum */
                    g++; /* Add one to g -- another to get rid of */
                } 
            }
        }
        h++;
        if(g == 0){  /* If cycles all through the for loop above, but */
            break; /* nothing changes, then no need to repeat */
        } 
    } while(h < Liv); /* Once the loop cycles through without changes, stops */
    /*Below brings back any individuals that are parents from -1 (removed) */
    for(i=0; i<Liv; i++){ 
        if(ID[i][4] < 0){ /* For all individuals in ID, find the dead ones */
            for(j=0; j<Liv; j++){ /* Check again to see if they are parents in ID */
                if(ID[i][0]==ID[j][5] || ID[i][0]==ID[j][6]){
                    ID[i][4] = M+1; /* If so, make very old */
                } /* Older than can do anything */
            } /* Effectively flagged as existing only to calculate kinship */
            for(j=0; j<l; j++){ /* Same for checking if parents in OFF */
                if(ID[i][0]==OFF[j][5] || ID[i][0]==OFF[j][6]){
                    ID[i][4] = M+1;
                }
            } /* Living dead (retained because parent of living) */
        } /* Living dead cannot mate, reproduce, etc. Can only be used in kin calcs */
        if(ID[i][4] < -1){ /* If became -2, make -1 for later removal */
            ID[i][4] = -1;
        } /* Now all removable individuals are -1 */
    }
}

