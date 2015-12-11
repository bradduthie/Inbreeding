/***************************************************************/
/* Program: randunif.c
   By: Brad Duthie                                             
   Description: Will sample size "samp" from Vec without replace
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***************************************************************/

#include "randunif.h" /* Need random unform numbers */
#include "array.h" /* Array used to handle vectors */

void SampleWithoutReplacement(int Pop, int samp, int *Vec){

    int t = 0; /* total input records dealt with */
    int m = 0; /* number of items selected so far */
    
    double u;

    while (m < samp){
        u = randunif();
        if ((Pop - t)*u >= samp - m ){
            t++;
        }
        else{
           Vec[m] = t;
           t++; 
           m++;
        }
    }
}
