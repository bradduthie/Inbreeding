#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void RegrTr(double **REGR, double **OFF, int prevOff);

void inbrdep(double **ID, double **REGR, int prevOff, int Liv, int M);

void printMeanKinship(double **Rmof, double **ID, int Liv, int rep, int gen, double Beta1,
    double cost, int DOM);
