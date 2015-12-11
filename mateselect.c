/***********************************************************************************/
/* Program: mateselect.c
   By: Brad Duthie                                             
   Description: Allows the dominant sex to select a mate based on a set of criteria.
   Compile: gcc mateselect.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void mateselect(double **ID, double **Rmat, int Liv, int DOM, int Nloci, int M, int mc,
    int Active, double alpha, int avoid, int Kind){
    /* Note that ID should be sorted now, dominant sex individuals at array top */    
    /* Note that dominant individuals should be sorted by quality at top of ID */

    int i, j, k, h, ch, s, count, no, sel, choo;
    double r, new, *PrM, *PrT, PrS, PrR, PrC;

    /* ===========================================================*/
    /* Determine where DOM and non-DOM sex row begins in ID array */
    /* ===========================================================*/        
    k = 0; /*Below counts the number of the dominant sex */
    for(i=0; i<Liv; i++){
        if(ID[i][1]==DOM && ID[i][4] > -1 && ID[i][4] <= M && ID[i][2] > -1){
            k++; 
        }else{
            break;
        }
    } /* So the dominants that can breed end at k */

    s = 0; /* Below gets the starting point for non-doms */
    for(i=0; i<Liv; i++){
        if(ID[i][1]==DOM){
            s++;
        }else{
            break;
        }
    } /* So the non-dominant sex starts at s */

    /* =======================================================================*/
    /* Male adjustment -- cannot allow to select more females than available  */
    /* =======================================================================*/
    if(DOM == 1){ /* If males are the choosing sex */
        count = 1; /* count is 1 because females can only be selected once */
        if(mc > (Liv-k)){ /* But the times males select needs to be determined */
            sel = Liv - k; /* This is no more than the total female number */
        }else{ 
            sel = mc; /* Or, more likely, just the parameter mc */
        }
    }else{
        count = mc; /* If females are choosing sex, males can be selected mc times */
    }

    /* =======================================================================*/
    /* Each DOM individual assesses and chooses one or more mates             */
    /* =======================================================================*/

    MAKE_1ARRAY(PrM,Liv); /* Vector for the probability of selecting an individual */
    MAKE_1ARRAY(PrT,Liv); /* Cumulative probability vector used for sampling */

    for(i=0; i<k; i++){ /* For each member of dominant sex, starting with best */
        for(j=0; j<Liv; j++){ /* For each other individual in the population */
            PrM[j] = 0; /* Just start everything off with zero */
            PrT[j] = 0; /* The cumulative vector elements are zeros too */
        }
        ch  = -1; /* Initially, the individual has not made a mate choice(s) */
        no  = 0;  /* Used to break out of loop later if no more mates available */
        for(j=s; j<Liv; j++){ /* Check out each of opposite sex available */
            /* If potential mate not taken too many times or dead */
            if(ID[j][1]>(-1*count) && ID[j][4]>-1 && 
            ID[j][4]<=M && ID[j][2] >= 0 && ID[j][3]>=0){
                if(Kind == 1){ /* If complete kin recognition */
                    r = Rmat[i][j+1]; /* Kinship with potential mate */
                }else{/* Below restricts to siblings */
                    if(ID[i][5]==ID[j][5] || ID[i][6]==ID[j][6]){
                        r = Rmat[i][j+1];
                    }else{ /* If not a sibling, treat as unrelated */
                        r = 0;
                    }
                }
                /* Below compares existing best mate (o), to new (o) */
                /* More complex rules using mate traits can be added */
                new = 1; /* Default quality of 1 -- adjusted with kinship */
                for(h=0; h<Active; h++){ /* Strategy alleles affect quality */
                    if(ID[i][(h*4)+10]==0){
                        new += r*alpha;
                    }
                    if(ID[i][(h*4)+11]==0){
                        new += r*alpha;
                    }
                } /* new is now higher if a relative */
                if(avoid == 1){ /* 1/new for avoidance, new for pref */
                    PrM[j] = 1/new; 
                }else{
                    PrM[j] = new;
                }
                no++; /* Adds recognition that a mate(s) was available */
            } /* Quality now increased (new) or decreased (1/new) with kinship */
        }
        if(no == 0){ /* Checks to see if a mate was available in that last loop */
            break; /* No need to continue if there are no mates left. */
        }
    
        /* ================================================================ */
        /* SELECTION IS PROBABILISTIC BASED ON MATE QUALITY                 */
        /* ================================================================ */
        switch(DOM){ /* Need to do different things depending on choosing sex */
            case 0: /* How to select if females are choosing */
                PrS = 0; /* Sum of the prob vector for choosing */
                for(j=0; j<Liv; j++){
                    PrS += PrM[j];
                } /* PrS now can be thought of as `total quality' */
                if(PrS > 0){ /* If some probability of getting mate */
                    PrR = 0; /* Running sum of probability vect */
                    for(j=0; j<Liv; j++){
                        PrT[j] = (PrM[j] / PrS) + PrR;
                        PrR = PrT[j];
                    } /* PrT now increasing vector, elements 0-1 */
                    PrC = randunif(); /* Vector prob select*/
                    j = 0; /* j increases to find position selected */
                    while(PrT[j] < PrC){
                        j++; /* Stops when random PrC pos found */
                    }
                    ch = j; /* The choice of individual is now j */
                    PrS -= PrM[j];
                    PrM[j] = 0;
                    ID[i][8] = ch; /* Add i's choice into her col 8 */
                }else{ /* Else there is no probability of getting mate */
                    ch = -1; /* So use -1 to indicate no mate acquired */
                }    
                if(ch > -1){ /* If i actually made some choice */
                    ID[ch][8] = i; /* Change the mate */
                    if(ID[ch][1] >= 0){ /* If mate not claimed  */
                        ID[ch][1] = -1; /* Mate 1st time */
                    }else{
                        ID[ch][1]--; /* Else subtract one */
                    } /* Above tracks times ch has been selected */
                    ID[ch][2] = ID[i][2]; /* non to dom */    
                    ID[ch][3] = ID[i][3]; /* non-dom to dom */
                }
                break; 
            case 1: /* Selection if males are choosing */
                PrS  = 0; /* Sum of the prob vector for choosing */
                choo = 0; /* choices a male can make */
                for(j=0; j<Liv; j++){
                    PrS += PrM[j];
                    if(PrM[j] > 0){
                        choo++; /* Finds number potential females */
                    }
                } /*If male can select more than avail, just mates with all */
                if(sel >= choo){
                    for(j=0; j<Liv; j++){
                        if(PrM[j] > 0){
                            ID[j][8] = i;
                            ID[j][1] = -1;
                            ID[j][2] = ID[i][2];
                            ID[j][3] = ID[i][3];
                            ID[i][8] = j;
                        }
                    }
                }else{ /* If their are more females than he can select */
                    choo = sel; /* Prob selection, as with females */
                    while(choo > 0){ /* But with choo choices, not 1 */
                        PrR = 0;/*Running sum of prob vector*/
                        for(j=0; j<Liv; j++){
                            PrT[j] = (PrM[j] / PrS) + PrR;
                            PrR = PrT[j];
                        }
                        PrC = randunif();
                        j = 0;        
                        while(PrT[j] < PrC){
                            j++;
                        }
                        ID[j][8] = i;
                        ID[j][1] = -1;
                        ID[j][2] = ID[i][2];
                        ID[j][3] = ID[i][3];
                        ID[i][8] = j;
                        choo--;
                    }
                }
                break;
            default: /* This just prints out if for some reason DOM != {0,1} */
                printf("%d\tERROR IN MATE SELECT\n",DOM);
                break;
        }
    } 
    FREE_1ARRAY(PrM); /* Free the memory; vector is no longer needed */
    FREE_1ARRAY(PrT); /* Free the memory; vector is no longer needed */
}


