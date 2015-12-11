/***********************************************************************************/
/* Program: arraysort2D.c
   By: Brad Duthie                                             
   Description: Sorts a 2D array by the value of a column element.
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include "array.h"

void vectorrank(double *vec, double **Rank, int xlen){
    
    int i, j;
    double *r, m, n;

    MAKE_1ARRAY(r,xlen*xlen);


    for(i=0; i<(xlen*xlen); i++){
        r[i] = 0;
        for(j=0; j<(xlen*xlen); j++){
            if(vec[i] > vec[j]){
                r[i]++;
            }    
        }
    }

    m = 0.0; /*row */
    n = 0.0; /* column */
    for(i=0; i<(xlen*xlen); i++){
        Rank[i][0] = m;
        Rank[i][1] = n;
        Rank[i][2] = r[i];
        m++;
        if(m==xlen){
            n++;
            m = 0;
        }
    }

    FREE_1ARRAY(r);
    
}


