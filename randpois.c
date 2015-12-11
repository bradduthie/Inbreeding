/***************************************************************/
/* Program: randpois.c
   By: Brad Duthie                                             
   Description: Will produce a Poisson distributed random number
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***************************************************************/


#include "randunif.h" /* Need randunif header file */

int randpois(double lambda){
    double L, u;
    double p = 1;
    int k = 0;
    L = exp(-1*lambda);
    while(p > L){
        k++;
        u = randunif();
        p = p * u;
    }
    return k-1;
}
