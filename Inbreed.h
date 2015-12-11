#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "array.h"
#include "randnormINT.h"
#include "randnorm.h"
#include "landscape.h"
#include "arraysort2D.h"
#include "vectorrank.h"
#include "Rbuild.h"
#include "mateselect.h"
#include "babies.h"
#include "parents.h"
#include "mortality.h"
#include "retain.h"
#include "deadbottom.h"

#define CORES 2
#define length(x) (sizeof (x) / sizeof *(x))

void Inbreed(int mc, int M, int Imm, int Clu, int DOM, double *RES, double Beta1, int rep,
    int Active, int Neutral, int load, double cost, double alpha, int avoid,
    int gen, int muSt, double mu, double hval, int Kind, int xlen, int mNalleles,
    int sdNalleles, int prP);
