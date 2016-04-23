#include <stdio.h>
#include "stdlib.h"
#include <complex.h>
#include <math.h> 

#define TRUE 1
#define FALSE 0

typedef struct _COMPLEX {
   double real, imag;
} COMPLEX;

int Powerof2(int n,int *m,int *twopm);
int FFT2D(COMPLEX **c,int nx,int ny,int dir);
int FFT(int dir,int m,double *x,double *y);

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
int FFT2D(COMPLEX **c,int nx,int ny,int dir) {
   int i,j;
   int m,twopm;
   double *real,*imag;

   /* Transform the rows */
   real = (double *)malloc(nx * sizeof(double));
   imag = (double *)malloc(nx * sizeof(double));
   if (real == NULL || imag == NULL)
      return(FALSE);
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return(FALSE);
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = c[i][j].real;
         imag[i] = c[i][j].imag;
      }
      FFT(dir,m,real,imag);
      for (i=0;i<nx;i++) {
         c[i][j].real = real[i];
         c[i][j].imag = imag[i];
      }
   }
   free(real);
   free(imag);

   /* Transform the columns */
   real = (double *)malloc(ny * sizeof(double));
   imag = (double *)malloc(ny * sizeof(double));
   if (real == NULL || imag == NULL)
      return(FALSE);
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return(FALSE);
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = c[i][j].real;
         imag[j] = c[i][j].imag;
      }
      FFT(dir,m,real,imag);
      for (j=0;j<ny;j++) {
         c[i][j].real = real[j];
         c[i][j].imag = imag[j];
      }
   }
   free(real);
   free(imag);

   return(TRUE);
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
int FFT(int dir,int m,double *x,double *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;
   double u[3], t[3], c[3]; 
   int kl = 1; 

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;
      //printf("Here\n"); 
   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
	 kl++; 	
         j -= k;
         k >>= 1;
      }
      j += k;
      printf("Here1\n"); 
      printf("%d\n", kl); 
   }
    printf("Here2\n"); 

   /* Compute the FFT */
   c[1] = -1.0;
   c[2] = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u[1] = 1.0;//create an array for these 
      u[2] = 0.0;
      u[1] = 0.0; 
      u[2] = 1.0; 
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t[1] = u[1] * x[i1] - u[2] * y[i1];
            t[2] = u[1] * y[i1] + u[2]  * x[i1];
            x[i1] = x[i] - t[1];
            y[i1] = y[i] - t[2];
            x[i] += t[1];
            y[i] += t[2];
         }
         z =  u[1] * c[1] - u[2] * c[2];
         u[2] = u[1] * c[2] + u[2] * c[1];
         u[1] = z;
      }
      c[2] = sqrt((1.0 - c[1]) / 2.0);
      if (dir == 1)
         c[2] = -c[2];
      c[1] = sqrt((1.0 + c[1]) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return(TRUE);
}

/*-------------------------------------------------------------------------
   Calculate the closest but lower power of two of a number
   twopm = 2**m <= n
   Return TRUE if 2**m == n
*/
int Powerof2(int n,int *m,int *twopm)
{
   if (n <= 1) {
      *m = 0;
      *twopm = 1;
      return(FALSE);
   }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
      return(FALSE);
   else
      return(TRUE);
}

int main(int argc, char *argv[]) {
	
	int n = 8; 	
	int poweroftwo = 8; 
	int m = 1; 
	printf("%d\n", Powerof2(n, &m, &poweroftwo)); 
	double test[n]; 
	double test1[n]; 
	for (n = 0; n < 8; n++){
		test[n] = 1.0; 
		test1[n] = 2.1; 
	} 
	
	FFT(1, 8, test, test); 
	printf("Gets here\n"); 

	
   
}
