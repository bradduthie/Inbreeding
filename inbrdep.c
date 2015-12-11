/***********************************************************************************/
/* Program: inbrdep.c 
   By: Brad Duthie  
   
   Description: Two functions to transfer offspring information and check survival 
   Compile: gcc inbrdep.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

/* This function transfers information from OFF to REGR, as needed */
void RegrTr(double **REGR, double **OFF, int prevOff){

    int i;

    for(i=0; i<prevOff; i++){
        REGR[i][0] = OFF[i][0]; /* Get individual number */
        REGR[i][1] = OFF[i][7]; /* This will store an individual f coefficient */
        REGR[i][2] = 0.5 * (2 - OFF[i][10] + OFF[i][11]); /* p-score for allele */
        REGR[i][3] = OFF[i][1];
        REGR[i][4] = 1; /* Start assuming complete survival probability, adj later */
        REGR[i][5] = 0; /* Assume dead unless proved otherwise by becoming an adult */
        REGR[i][6] = OFF[i][2]; /* This was the juvenile probability of survival */
        REGR[i][7] = 0; /* An individual's offspring production will go here */
    }
}

/* This function checks to see if individuals from original OFF ultimately survived */
void inbrdep(double **ID, double **REGR, int prevOff, int Liv, int M){

    int i, j;

    for(i=0; i<Liv; i++){
        if(ID[i][4] > -1 && ID[i][4] <= M){ /* If i is alive */
            for(j=0; j<prevOff; j++){ /* Go through REGR and find it, j */
                if(ID[i][0] == REGR[j][0]){
                    REGR[j][5] = 1; /* Add a one to indicate survival */
                    break; /* No need to keep through this loop */
                }
            }
        }
    }
}


/* Prints kinship relative to an indivual to a kinship file (use with pedigree) */
void printMeanKinship(double **Rmof, double **ID, int Liv, int rep, int gen, double Beta1,
    double cost, int DOM){

    int    l, j, k, h, subs, t, m, *ksamp, l90p, l99p, u90p, u99p;
    double meanK, /*meanI,*/ SDk, *kvals, LQ, UQ, medk, lowk, highk, U90, L90, U99, L99;
    double dtemp, u, *ksub;
    /* FILE *kinshipmeans; */   /* This will print out if Pedigree is also printing */
    FILE *kinshipsraw;    /* Will print out raw kinships*/
    /*FILE *kinshipssampl;*/  /* Will print out random sample kinships*/


    /* kinshipmeans = fopen("kinshipmeans.txt","a+"); */
    kinshipsraw   = fopen("kinshipsraw.txt", "a+");
    /*kinshipssampl = fopen("kinshipsampl.txt", "a+");*/

    meanK = 0.0;
    l     = 0;
    subs  = 10000;

    fprintf(kinshipsraw,"%d\t%f\t%f\t%d\t",DOM,cost,Beta1,gen);
    for(j=0; j<Liv; j++){ /* This loop gets the mean k between individuals */
        for(k=0; k<Liv; k++){
            if(j != k && ID[j][4]==0 && ID[k][4]==0){ /* (k+1) */
                meanK += Rmof[j][k+1];    
                l++;
            }
        }
    }
    meanK = meanK / l;

    if(subs > l){
        subs = l;
    }

    MAKE_1ARRAY(kvals,l);
    MAKE_1ARRAY(ksamp,subs);
    MAKE_1ARRAY(ksub,subs);

    h = 0;
    for(j=0; j<Liv; j++){ /* This loop gets the mean k between individuals */
        for(k=0; k<Liv; k++){
            if(j != k && ID[j][4]==0 && ID[k][4]==0){ /* (k+1) */
                SDk      += (Rmof[j][k+1]-meanK)*(Rmof[j][k+1]-meanK);
                kvals[h] =  Rmof[j][k+1];
                h++;
            }
        }
    }
    SDk = sqrt(SDk / l);

    t     = 0;
    m     = 0;

    while (m < subs){
        u = randunif();
        if ((l - t)*u >= subs - m ){
            t++;
        }
        else{
           ksamp[m] = t;
           t++; 
           m++;
        }
    }

    /*fprintf(kinshipssampl,"%d\t",gen); */
    for(j=0; j<subs; j++){
        ksub[j] = kvals[ksamp[j]];
        /*fprintf(kinshipssampl,"%f\t",ksub[j]);*/
    }
    /*fprintf(kinshipssampl,"\n");*/


    for(j=subs-2; j>=0; j--){ 
        for(k=0; k<=j; k++){ 
            if(ksub[k] > ksub[k+1]){ 
                dtemp     = ksub[k];
                ksub[k]   = ksub[k+1];
                ksub[k+1] = dtemp;
            }
        }     
    }

    l90p = floor(subs*0.1);
    u90p = floor(subs*0.9);
    l99p = floor(subs*0.01);
    u99p = floor(subs*0.99);

    lowk  = ksub[0];
    highk = ksub[subs-1];
    medk  = ksub[subs/2];
    LQ    = ksub[subs/4];
    UQ    = ksub[(subs/2) + (subs/4)];
    U90   = ksub[u90p];
    L90   = ksub[l90p];
    U99   = ksub[u99p];
    L99   = ksub[l99p];

    fprintf(kinshipsraw,
            "%d\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            DOM,cost,Beta1,gen,l,meanK,SDk,lowk,highk,medk,LQ,UQ,L90,U90,L99,U99);
    
    /* This loop prints relevant information regarding kinship */
    /*for(j=0; j<Liv; j++){
        if(ID[j][4]==0){ 
            fprintf(kinshipmeans,"%d\t%d\t%f\t",rep,gen,ID[j][0]);
            fprintf(kinshipmeans,"%f\t",ID[j][1]); 
            fprintf(kinshipmeans,"%f\t",ID[j][5]); 
            fprintf(kinshipmeans,"%f\t",ID[j][6]); 
            fprintf(kinshipmeans,"%f\t",(ID[j][10]+ID[j][11]));
            meanI = 0.0;
            l     = 0;
            for(k=0; k<Liv; k++){
                if(j != (k+1) && ID[k][4]==0){
                    meanI += Rmof[j][k+1];
                    l++;
                }
            } 
            fprintf(kinshipmeans,"%f\t%f\t%f\n",meanI/l,meanK,SDk);
        }
    } */ /* Prints an individual's mean kinship with others, overall mean in pop */

    FREE_1ARRAY(kvals);
    FREE_1ARRAY(ksamp);
    FREE_1ARRAY(ksub);
    /*fclose(kinshipmeans);*/
    fclose(kinshipsraw);
    /*fclose(kinshipssampl);*/
}








