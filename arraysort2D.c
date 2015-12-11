/***********************************************************************************/
/* Program: arraysort2D.c
   By: Brad Duthie                                             
   Description: Sorts a 2D array by the value of a column element.
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include<stdio.h>

void arraysort2D(double **ID, int rows, int cols, int by, int incr){
    
    int j, k, h;
    double dtemp;

    switch (incr) {
        /*-------------------------------*/
        case 1: /* If 'incr' is 1 (true) */
        /*-------------------------------*/
        for(j=rows-2; j>=0; j--){ 
            for(k=0; k<=j; k++){ 
                if(ID[k][by] > ID[k+1][by]){
                    for(h=0; h<cols; h++){
                        dtemp = ID[k][h];
                        ID[k][h] = ID[k+1][h];
                        ID[k+1][h] = dtemp;
                    }
                }     
            }
        }                        
        break;
        /*-------------------------------*/
        default: /* If 'incr' is not 1 (false) */
        /*-------------------------------*/
        for(j=rows-2; j>=0; j--){ 
            for(k=0; k<=j; k++){ 
                if(ID[k][by] < ID[k+1][by]){
                    for(h=0; h<cols; h++){
                        dtemp = ID[k][h];
                        ID[k][h] = ID[k+1][h];
                        ID[k+1][h] = dtemp;
                    }
                }     
            }
        }
        break;
    }
}


