#include <time.h>
#include <complex.h>
#include "FFT.h"
#include "randunif.h"
#include "array.h"


#define Pi 3.14159265358979323846264338327


/* Define complex numbers */
typedef struct _complex{  
       double real;  
       double img;
} Complex;
/* Can now use complex numbers */

Complex PrepFFT(double S_f, double phi){  
   Complex x2;  
   x2.real = pow(S_f,0.5) * (cos(2.0*Pi*phi));  
   x2.img  = pow(S_f,0.5) * (sin(2.0*Pi*phi));  
   return x2;
 }  

int landscape(double *LS, double BETA, int X){

    int i, j;
    int LG = log(X*X)/log(2);
    
    size_t x = X;
    size_t L = X*X;
    
    double *phi, **u, **v, **S_f, **p, *w, *y;
    double z, s;
    
    Complex **q; 

    MAKE_1ARRAY(phi,L);
    MAKE_1ARRAY(w,L);
    MAKE_1ARRAY(y,L);
    MAKE_2ARRAY(u,x,x);
    MAKE_2ARRAY(v,x,x);
    MAKE_2ARRAY(S_f,x,x);
    MAKE_2ARRAY(p,x,x);
    MAKE_2ARRAY(q,x,x); 
    
    /* Below makes u and its trandpose, v */
    for(i=0; i<X; i++){
        z = 0;
        for(j=0; j<(X/2); j++){
            u[i][j] = z/X;
            v[j][i] = z/X;
            z++; 
        }
        for(j=(X/2); j<X; j++){
            u[i][j] = z/X;
            v[j][i] = z/X;
            z--; 
        }        
    } /* Now have matrices u and v */
    
    
    /* Get the matrix S_f */
    for(i=0; i<X; i++){
        for(j=0; j<X; j++){
            s = (u[i][j]*u[i][j] + v[i][j]*v[i][j]);
            if(s > 0){
                S_f[i][j] = pow(s,(BETA/2));
            }else{
                S_f[i][j] = 0;    
            }
        }    
    } /*Matrix S_f made from u and v */


    /*Below makes random unifs, phi */
    for(i=0; i<L; i++){
        phi[i] = randunif();
    }
    
    /* Puts the random unifs in a matrix */
    for(i=0; i<X; i++){
        for(j=0; j<X; j++){
            p[i][j] = phi[(X*i)+j];
            q[i][j] = PrepFFT(S_f[i][j],p[i][j]);
            w[(X*i)+j] = q[i][j].real;
            y[(X*i)+j] = q[i][j].img;            
        }
    } /* Matrix made */
    
    /*Below runs the fast Fourier transformation */
    FFT(1,LG,w,y);

    /*Below puts this in the vector LS */
    for(i=0; i<X; i++){
        for(j=0; j<X; j++){
            LS[(X*i)+j] = w[(X*i)+j];
        }
    }
    
    FREE_1ARRAY(w);
    FREE_1ARRAY(y);
    FREE_1ARRAY(phi);
    FREE_2ARRAY(u);    
    FREE_2ARRAY(v);
    FREE_2ARRAY(S_f);    
    FREE_2ARRAY(p);    
    FREE_2ARRAY(q);                    
    
    return EXIT_SUCCESS;    
}

