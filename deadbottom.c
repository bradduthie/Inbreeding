/***********************************************************************************/
/* Program: deadbottom.c
   By: Brad Duthie                                             
   Description: Brings positive values to matrix top
   Compile: gcc deadbottom.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include<stdio.h>

void deadbottom(double **ID, int rows, int cols, int M){
    
    int j, k, h;
    double dtemp;
    
    /**********************************************************************/
    /* Below shifts all the dead, and effectively dead (too old to breed) */
    /* individuals from the population to the bottom of the table */
    /**********************************************************************/
    for(j=rows-2; j>=0; j--){ 
        for(k=0; k<=j; k++){ 
            if( (ID[k][4]<0 || ID[k][4]>M) && 
                (ID[k+1][4]>=0 && ID[k+1][4] <= M) ){
                for(h=0; h<cols; h++){
                    dtemp = ID[k][h];
                    ID[k][h] = ID[k+1][h];
                    ID[k+1][h] = dtemp;
                }
            }     
        }
    }
}


