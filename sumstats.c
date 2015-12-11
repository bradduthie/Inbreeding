/***********************************************************************************/
/* Program: sumstats.c
   By: Brad Duthie                                             
   Description: Returns summary stats for the inbreed program.
   Compile: gcc mortality.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"

void sumstats(double *RES, double **OFF, int M, int Nloci, int loadstart, int load, 
    int Neutstart, int Neutral, int Active, int l, int Imm, int xlen, int DOM,
    double **REGR, int generation, int prevOff){

    int i, j;    
    double g, k, h, m, b, dels;
    double sumxy, sumx, sumy, sumx2;    


    /* ==========================================================*/
    /* Calculate the 0 allele frequency =========================*/
    /* ==========================================================*/        
    g    = 0;
    k    = 0;
    h    = 0;
    dels = 0;
    for(i=0; i<l; i++){
        if(OFF[i][4] >= 0 && OFF[i][4] <= M){
            for(j=0; j<Active; j++){
                if(OFF[i][((4*j)+10)] == 0){
                    g++;
                }
                if(OFF[i][((4*j)+11)] == 0){
                    g++;
                }
            }
            for(j=0; j<Neutral; j++){
                if(OFF[i][((4*j)+Neutstart)] == 0){
                    h++;
                }
                if(OFF[i][((4*j)+Neutstart+1)] == 0){
                    h++;
                }
            }
            for(j=0; j<load; j++){
                if(OFF[i][((4*j)+loadstart)] == 0){
                    dels++;
                }
                if(OFF[i][((4*j)+loadstart+1)] == 0){
                    dels++;
                }
            }
            k += 2;
        }
    }

    RES[0] = g / (Active * k); /* Strategy allele frequency */
    RES[1] = h / (Neutral * k); /* Neutral allele frequency */
    RES[2] = dels / (load * k); /* Load allele frequency */

    /* ==========================================================*/
    /* Regression coefficients -- juvenile survival =============*/
    /* ==========================================================*/    

    if(generation > 0){
        sumx   = 0; /* Sum of inbreeding coefficients */
        sumy   = 0; /* Sum of probability of survival */
        sumxy  = 0; /* Sum of inbreeding coefficient times prob survival */

        for(i=0; i<prevOff; i++){
            sumx   += REGR[i][1]; /* Individual's f coefficient */
            sumy   += log(REGR[i][6]); /* LOG probability of juvenile survival */
            sumxy  += REGR[i][1] * log(REGR[i][6]);
            sumx2  += REGR[i][1] * REGR[i][1];
        }

        /* Computing the slope below */
        m = (prevOff * sumxy - sumx * sumy) / (prevOff * sumx2 - sumx * sumx);
        /* Computing the intercept below */
        b = (sumy * sumx2 - sumx * sumxy) / (prevOff * sumx2 - sumx * sumx);
    
        RES[3]  = sumx / prevOff; /* Mean inbreeding coefficient */
        RES[4]  = sumy / prevOff; /* Mean log probability of juvenile survival */
        RES[5]  = m; /* Slope of juvenile inbreeding depression regression */
        RES[6]  = b; /* Intercept of juvenile inbreeding depression regression */
    }

    /* ==========================================================*/
    /* An inbreeding depression (sort of) estimate using survival*/
    /* ==========================================================*/    

    if(generation > 0){
        sumx   = 0; /* Sum of inbreeding coefficients */
        sumy   = 0; /* Sum of probability of survival */
        sumxy  = 0; /* Sum of inbreeding coefficient times prob survival */

        for(i=0; i<prevOff; i++){ 
            sumx   += REGR[i][1]; /* What was the individual's f coefficient? */
            sumy   += REGR[i][5]; /* Did the indivual survive (1) or not (0) */
            sumxy  += REGR[i][1] * REGR[i][5];
            sumx2  += REGR[i][1] * REGR[i][1];
        }

        /* Computing the slope below */
        m = (prevOff * sumxy - sumx * sumy) / (prevOff * sumx2 - sumx * sumx);
        /* Computing the intercept below */
        b = (sumy * sumx2 - sumx * sumxy) / (prevOff * sumx2 - sumx * sumx);
    
        RES[7]  = m; /* This will later be lifetime inbreeding depression slope */
        RES[8]  = b; /* This will later be lifetime inbreeding depression intercept */
    }

    /* ==========================================================*/
    /* An inbreeding depression estimate with reproductive output*/
    /* ==========================================================*/    

    if(generation > 0){
        /* All of this will calculate the old offspring's ultimate reproduction */
        for(i=0; i<prevOff; i++){ /* Go through each of the offspring from last gen */
            if(REGR[i][5] == 1){ /* If they were alive, might have had offspring */
                for(j=0; j<l; j++){ /* So go through OFF and check */
                    if(REGR[i][0]==OFF[j][5] || REGR[i][0]==OFF[j][6]){
                        REGR[i][7]++; /* Add existing offspring */
                    }
                }
            }
        } /* Now REGR[i][7] should have the offspring of each individual */

        /* We use these offspring production calculations for another ID measure */
        sumx   = 0; /* Sum of inbreeding coefficients */
        sumy   = 0; /* Sum of offspring count */
        sumxy  = 0; /* Sum of inbreeding coefficient times offspring count */

        for(i=0; i<prevOff; i++){ 
            sumx   += REGR[i][1]; /* What was the individual's f coefficient? */
            sumy   += REGR[i][7]; /* How many offspring did the individual have? */
            sumxy  += REGR[i][1] * REGR[i][7];
            sumx2  += REGR[i][1] * REGR[i][1];
        }

        /* Computing the slope below */
        m = (prevOff * sumxy - sumx * sumy) / (prevOff * sumx2 - sumx * sumx);
        /* Computing the intercept below */
        b = (sumy * sumx2 - sumx * sumxy) / (prevOff * sumx2 - sumx * sumx);
    
        RES[9]   = m; /* This will later be lifetime inbreeding depression slope */
        RES[10]  = b; /* Will later be lifetime inbreeding depression intercept */
    }

    /* ==========================================================*/
    /* Write over RES[9] and RES[10] with load sd, population    */
    /* ==========================================================*/    

    k = 0;
    b = 0;
    for(i=0; i<l; i++){
        if(OFF[i][4] >= 0 && OFF[i][4] <= M){
            for(j=0; j<load; j++){
                if(OFF[i][((4*j)+loadstart)] == 0){
                    dels++;
                }
                if(OFF[i][((4*j)+loadstart+1)] == 0){
                    dels++;
                }
            }
            m =  dels / (2 * load);
            b += (m-RES[2])*(m-RES[2]);
            k++;
        }
    }
    b = sqrt((1/k)*b);
    
    RES[9]  = b;
    RES[10] = k;
}





















