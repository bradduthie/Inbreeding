#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void mortality(double **ID, double **OFF, int Nloci, int Liv, int l, int M, 
    double **SCAPE, double Beta1, int DOM, int loadstart, int load, double cost,
    double hval);
