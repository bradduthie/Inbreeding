/***********************************************************************************/
/* Program: babies.c
   By: Brad Duthie                                             
   Description: Uses the tables ID and Roff to make new offspring.
   Compile: gcc babies.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"

int babies(double **ID, int Clu, int i, int DOM, int Liv){

    int Off,  j, tOff;

    Off  = 0;

    switch(DOM){
        /*===========================================================*/
        /* Calculates offspring produced for a nest if fem dom ======*/
        /*===========================================================*/
        case 0:
            if(ID[i][8]<0){    /* If female did not find a mate  */
                Off = 0;       /* She produces no offspring      */
            }else{             /* If she did find a mate         */
                Off = Clu;     /* She produces Clu offspring     */
            }
            break;
        /*===========================================================*/
        /* Calculates offspring produced for a nest if male dom =====*/
        /*===========================================================*/
        case 1:
            for(j=0; j<Liv; j++){
                tOff = 0;
                if(ID[j][8]==i){      /* Checks how many females male found */
                    tOff += Clu;      /* Adds Clu for each female he had */
                    ID[j][9] = tOff;
                }
                Off += tOff;          /* Off is the total number of -- Clu * Females */
            }
            break;
        default:
            printf("ERROR IN BABIES.C\n");
            break;
    }
    return Off;
}


