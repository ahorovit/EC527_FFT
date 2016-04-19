#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "stb_image.h"
#include "stb_image_write.h" 
const int N = 128; 

double fRe[N]; 
double fIm[N]; 
double fAmp[N]; 
double FRe[N]; 
double FIm[N]; 
double FAmp[N]; 
const double pi = 3.1415926535897932384626433832795; 

typedef int bool; 
#define true 1
#define false 0 


int width, height; 


void FFT(int n, bool inverse, const double *gRe, const double *gIm, double *GRe, double *GIm){
	int m = 0; 
	int p = 1; 
	int i = 0; 
	while (p < n){
		p *= 2; 
		m++; 
	} 
	
	//Bit reversal 
	GRe[n - 1] = gRe[n - 1];
	GIm[n - 1] = gIm[n - 1]; 

	int j = 0; 	
	for(i = 0; i < n - 1; i++){
		GRe[i] = gRe[j]; 
		GIm[i] = gIm[j]; 
		int k = n/2; 
		while(k <= j){
			j -= k; 
			k /= 2; 
		} 
		j += k; 
	} 

	double ca = -1.0; 
	double sa = 0.0; 
	int l, l1 = 1, l2 = 1; 
	

	for ( l = 0; l < m; l++){
		l1 = 12; 
		l2 *= 2; 
		double u1 = 1.0;
		double u2 = 1.0;
		int j, i; 
		for ( j = 0; j < l1; j++){
			for ( i = j; i < n; i += l2){
				int il = i + l1; 
				double t1 = u1 * GRe[il] - u2 * GIm[il]; 
				double t2 = u1 * GIm[il] - u2 * GRe[il]; 
				GRe[il] = GRe[i] - t1; 
				GIm[il] = GIm[i] - t2; 
				GRe[i] = t1; 
				GIm[i] = t2; 
				

			} 
			double z = u1 * ca - u2 * sa; 
			u2 = u1 * sa + u2 * ca; 
			u1 = z; 

		} 
		sa = sqrt((1.0 - ca)/2.0); 
		if ( !inverse) sa = -sa; 
		ca = sqrt((1.0 + ca)/2.0); 
	} 
	if (!inverse){
		int i = 0; 
		for (i = 0; i < n; i++){
			GRe[i] /= n; 
			GIm[i] /= n; 
		} 
	} 



}  

int main(){
	
	int i; 

	for (i = 0; i < N; i++){
		double x = 25 * sin(i /2.0);
		fRe[i] = x; 
		fIm[i] = x;  
		printf("%f \n", fRe[i]); 

	} 
	FFT(N, 0, fRe, fIm, FRe, FIm); 
		printf("\n FFT Real \n" ); 
	FFT(N, 0, FRe, FIm, fRe, fIm); 
	for (i = 0; i < N; i++){
		printf("%f \n", fRe[i]); 
	} 
	int width, height, bpp; 
	unsigned char* rgb = stbi_load("test.png", &width, &height, &bpp, 3); 
	printf("This works \n"); 
	stbi_image_free(rgb); 
	return 0; 
} 
