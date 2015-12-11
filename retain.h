/***********************************************************************************/
/* Program: mortality.c
   By: Brad Duthie							         		
   Description: Marks individuals who are dead, but also parents (needed for Rmat)
   Compile: gcc mortality.c -ansi -Wall -pedantic	*/
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void retain(double **ID, double **OFF, int Liv, int l, int L);
