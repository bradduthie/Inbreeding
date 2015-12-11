/***************************************************************/
/* Program: randunif.c
   By: Brad Duthie                                             
   Description: Will produce a uniform distributed random number
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***************************************************************/

#include "as183.h"
/* as183() does the heavy lifting so "randunif()" can generate
   a single random uniform number */

double randunif(void){
    int seed[3];
    double randun;
       seed[0] = rand();
       seed[1] = rand();
       seed[2] = rand();
       randun  = as183(seed); 
    return randun;
}
