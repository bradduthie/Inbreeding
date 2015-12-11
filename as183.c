/***************************************************************/
/* Program: as183.c
   By: Brad Duthie							         		
   Description: Takes three seeds and returns a uniform random
   Compile: gcc as183.c -ansi -Wall -pedantic	*/
/***************************************************************/

#include<stdlib.h> /* Include standard library */
/* Below function returns the mod of abs(x)+1*/
/* This is used to verify an acceptable seed*/

int verify_seed(int x){
  x=abs(x) % 30000;
  return(++x);
 } /* Easy way of getting seeds */

double as183(int seeds[]){
  double unidev; /* Code below verifies the 3 seeds */
  seeds[0] = verify_seed(seeds[0]);
  seeds[1] = verify_seed(seeds[1]);
  seeds[2] = verify_seed(seeds[2]);
  /* Code below gets a decimal to be added to unidev */
  seeds[0] = (171 * seeds[0]) % 30269;
  seeds[1] = (172 * seeds[1]) % 30307;
  seeds[2] = (170 * seeds[2]) % 30323;
  /* unidev gets a random uniform number between zero and one */
  unidev = seeds[0]/30269.0 + seeds[1]/30307.0 + seeds[2]/30323.0;
  /* Return just the decimal, subtract integer of unidev */
  return(unidev - (int)unidev);
} /* We now have one random uniform number */
