/***********************************************************************************/
/* Program: InMult.c
   By: Brad Duthie                                             
   Description: Calls Inbreed.c multiple times to simulate different inbreeding
        scenarios.
   Compile: To compile with makefile, type `make' within directory, then hit ENTER
   Run:     To run, type `./InMult' on command line, then hit ENTER   */
/***********************************************************************************/

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Inbreed.h"
#include "array.h"

int main(void){

    int    mc, M, Imm, Clu, DOM, rep, i, j, avoid, xlen, Pedi, prP;
    int    Active, Neutral, load, gen, muSt, Kind, mNalleles, sdNalleles;
    double Beta1, cost, alpha, mu, hval, *RES, Beta1ch;

    /* =========== VARIABLES BETWEEN THE Xs BELOW ADJUST MODEL PARAMETERS ================*/
    /*  XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX   */
    /* ===================================================================================*/
    /* Model parameter values                                                             */
    /* ===================================================================================*/
    mc     = 5;      /* Number of females a male can mate with                            */
    M      = 0;      /* The age at which all individuals must die (must only use 0)       */
    Imm    = 5;      /* Number of immigrants per generation                               */
    Clu    = 5;      /* Clutch size for each female                                       */
    DOM    = 0;      /* Dominant sex (0 = female, 1 = male)                               */
    alpha  = 10.0;   /* Strength of inbreeding avoidance/preference                       */
    avoid  = 1;      /* Avoid (avoid = 1) or prefer (avoid = 0) allele mutant             */
    Beta1  = 0.0000; /* Selection coefficient on deleterious recessives                   */
    cost   = 0.0000; /* Cost of having an active inbreeding strategy                      */
    gen    = 3000;   /* Number of generations per replicate (default 3000)                */
    muSt   = 100;    /* Generation at which mutations may start                           */
    mu     = 0.001;  /* Mutation rate of any given allele                                 */
    hval   = 0.00;   /* The degree of partial dominance (0 = wild type completely dom)    */
    Kind   = 1;      /* Kin recognition (1 = recognise all; 0 = recognise only siblings)  */
    xlen   = 15;     /* Spatial x and y dimension -- carrying capacity = xlen*xlen        */
    /*                                                                                    */
    /* ===================================================================================*/
    /* Genome attributes of individuals                                                   */
    /* ===================================================================================*/
    Active     = 1;    /* Number of active alleles (preference or avoidance)              */
    Neutral    = 1;    /* Number of neutral alleles (no effect on fitness )               */
    load       = 1000; /* Number of inbreeding load alleles                               */
    mNalleles  = 2;    /* Set mean number of alleles per locus (must be >1)               */
    sdNalleles = 0;    /* Set SD of allele number per locus                               */
    /*                                                                                    */
    /* ===================================================================================*/
    /* Simulation details                                                                 */
    /* ===================================================================================*/
    rep        = 1;    /* Simulations run                                                */
    Pedi       = 0;     /* Print last rep's pedigree? (0:no, 1:yes) WARNING: 200Mb file   */
    Beta1ch    = 0.000; /* How much should Beta1 change each rep?                         */ 
    /*                                                                                    */
    /* ===================================================================================*/
    /*  XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX   */
    /* ===================================================================================*/

    srand(time(NULL) ^ (getpid()<<16)); /* Use time to generate random seed */

    i   = 0; /* Loops through different replicate simulations */
	prP = 0; /* Indicator for actually printing the pedigree */
    while(i < rep){
		if(i == rep-1 && Pedi == 1){ /* If it's the last rep, and should print pedigree */
			prP = 1;                 /* Switch this variable so pedigree will print */
		}
        prP = 2; /* Never printing XXX For now */
        MAKE_1ARRAY(RES,11); /* The RES array holds summary statistics */
        for(j=0; j<11; j++){ /* Need to refresh RES with zeros for each simulation */
            RES[j]=0.0;
        }

        /* The function below is the main simulation function from Inbreed.c */
        Inbreed(mc,M,Imm,Clu,DOM,RES,Beta1,i,Active,Neutral,load,cost,alpha,
            avoid,gen,muSt,mu,hval,Kind,xlen,mNalleles,sdNalleles,prP);
        /* Below prints all of the summary statistics into the terminal */
        printf("%f\t%d\t%f\t",Beta1,DOM,cost);
        for(j=0; j<11; j++){
            printf("%f\t",RES[j]);
        }
        printf("\n");
        FREE_1ARRAY(RES); /* Free the RES array after printing is finished */

        /* Below increases the selection coefficients for the next loop */
        Beta1 += Beta1ch;

        i++;

        /*if(i == 6 || i == 30){
            Beta1 = 0.0;
            cost  = 0.0025;
        }
        if(i == 12 || i == 36){
            Beta1 = 0.0;
            cost  = 0.005;
        }
        if(i == 18 || i == 42){
            Beta1 = 0.0;
            cost  = 0.01;
        }
        if(i == 24){
            DOM   = 1;
            Beta1 = 0.0;
            cost  = 0.0;
        } */

    }
    return 0;
}





