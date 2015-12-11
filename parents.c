/***********************************************************************************/
/* Program: parents.c
   By: Brad Duthie                                             
   Description: Determines which offspring came from which parents, and properties.
   Compile: gcc babies.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void parents(double **ID, double **OFF, double **Rmof, int *O, int Nloci, int Ind, 
    int k, int DOM, int l, int Clu, int Liv, int gener, int mc, double Beta1,
    double *RES, int rep, int loadstart, int load, int Active, int prP){

    int i, j, h, m, s, MateR, MomP, DadP, a1, a2;
    double deletes, iload, pscore, nscore; 
    FILE *Pedigree;

    /*===========================================================================*/
    /* First add offspring, assuming nestmates are parents ======================*/
    /*===========================================================================*/
    h = 0;
    for(i=0; i<l; i++){
        if(i==O[h]){ /* If offspring # hits count of off nest*/    
            h++; /* Move to the next nest */ 
        } 
        if(DOM == 0){ /* If females choosing sex */
            m = ID[h][8]; /* Row of Dad */ 
        }else{
            m = 0; /* Else males choosing sex, finds mums of nest h */
            while(h != ID[m][8]){ /* While m is not a mum in dad's nest */
                m++; /* Keep searching row of females */
            } /* When search is over, m is a mum of the nest */
            ID[m][9]--; /* One offspring of Mum used */
            if(ID[m][9] == 0){ /* Mum is done with male if clutch used */
                ID[m][8] = -1; /* Lose so don't find in loop 2nd time */
            } /* Now have found the first/next mum of the nest */
        } /* Or the male m that sired offspring in the nest */
        MateR     = ID[m][0];             /* Number of mate */
        OFF[i][0] = Ind;                  /* Label for the new individual born */
        Ind++;                            /* Add a new one to Ind */
        OFF[i][1] = floor(2*randunif());  /* Offspring randomly m or f */
        OFF[i][2] = ID[h][2];             /* Mother's nest coordinate x */
        OFF[i][3] = ID[h][3];             /* Mother's nest coordinate y */
        OFF[i][4] = 0;                    /* Starts at age zero */ 
        if(DOM == 0){
            OFF[i][5] = ID[h][0];     /* Offspring Mother Number */
            OFF[i][6] = MateR;        /* Offspring Father Number */
        }else{
            OFF[i][5] = MateR;        /* Mother Number */
            OFF[i][6] = ID[h][0];     /* Father Number */
        } /* Below gets Mum & Dad on Rmof matrix to calc f */
        MomP = 0; 
        DadP = 0; /* Below finds the matrix pos of Mom, Dad */
        while(OFF[i][5]!=Rmof[MomP][0]){
            MomP++; /* Keep increasing until MomP finds the Mum's position */
        }
        while(OFF[i][6]!=Rmof[DadP][0]){
            DadP++; /* Keep increasing until MomP finds the Mum's position */
        } /* Now have both Mum and Dad position in the kinship matrix `Rmof' */
        OFF[i][7] = Rmof[MomP][DadP+1]; /* Offspring inbreeding coefficient */
        OFF[i][8] = -1;  /* Offspring Mate selection -- no mate selected */
        OFF[i][9] = Clu; /* Offspring Clutch size */
        for(j=0; j<Nloci; j++){ 
            s = floor(randunif()*8); /* Randomly selects a number from 1 to 8 */
            /* There are eight possible ways to pass alleles parent-offspring */
            switch(s){ 
                case 0:
                    OFF[i][10+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][12+(j*4)];
                    break;
                case 1:
                    OFF[i][10+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][13+(j*4)];
                    break;
                case 2:
                    OFF[i][10+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][12+(j*4)];
                    break;
                case 3:
                    OFF[i][10+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][13+(j*4)];
                    break;
                case 4:
                    OFF[i][10+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][12+(j*4)];
                    break;
                case 5:
                    OFF[i][10+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][13+(j*4)];
                case 6:
                    OFF[i][10+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][12+(j*4)];
                    break;
                case 7:
                    OFF[i][10+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][13+(j*4)];
                    break;
                default:
                    printf("ERROR IN PARENTS.C\n");
                    break;
            }
        } /* 1 allele of each parent now randomly inherited by offspring at all loci */
    }

    /*===========================================================================*/
    /* PRINTS OFF PEDIGREE BELOW     ============================================*/
    /*===========================================================================*/

    if(prP == 1){
        Pedigree = fopen("Pedigree.txt","a+");
        for(i=0; i<l; i++){
            pscore = 0.0;
            nscore = 0.0;
            for(j=0; j<Active; j++){
                if(OFF[i][(j*4)+10]==0){
                    pscore++;
                }
                if(OFF[i][(j*4)+11]==0){
                    pscore++;
                }
            }
            for(j=0; j<Active; j++){
                if(OFF[i][(j*4)+14]==0){
                    nscore++;
                }
                if(OFF[i][(j*4)+15]==0){
                    nscore++;
                }
            }
            deletes = 0.0;
            iload   = 0.0;
            for(j=0; j<load; j++){
                a1 = OFF[i][((4*j)+loadstart)];
                a2 = OFF[i][((4*j)+loadstart)+1];
                if(a1 == 0){
                    iload++;
                }
                if(a1 == 0){
                    iload++;
                }            
                if(a1 == 0 && a2 == 0){
                    deletes++;
                }
            }

            pscore = pscore*0.5 / Active;
            nscore = nscore*0.5 / Active;

            fprintf(Pedigree,"%d\t%d\t%f\t",rep,gener,OFF[i][0]);
            fprintf(Pedigree,"%f\t%f\t%f\t",OFF[i][1],OFF[i][5],OFF[i][6]);
            fprintf(Pedigree,"%f\t%f\t",OFF[i][7],pscore);
            fprintf(Pedigree,"%f\t%f\t%f\t",RES[1],deletes,iload);
            fprintf(Pedigree,"%f\n",nscore);
        }
        fclose(Pedigree);
    }

}








