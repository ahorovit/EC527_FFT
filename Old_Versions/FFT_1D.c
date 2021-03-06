#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h" 
const int N = 128; 

double fRe[N][N][3]; 
double fIm[N][N][3]; 
double fAmp[N][N][3]; 
double FRe[N][N][3]; 
double FIm[N][N][3]; 
double FAmp[N][N][3]; 
const double pi = 3.1415926535897932384626433832795; 

typedef int bool; 
#define true 1
#define false 0 


int width, height, bpp; 


void FFT(int n, bool inverse, const double *gRe, const double *gIm, double *GRe, double *GIm, int stride, double factor){
	int m = 0; 
	int p = 1; 

	int i = 0; 

	while (p < n){
		p *= 2; 
		m++; 
	} 
	
	//Bit reversal 
	GRe[(n - 1) * stride] = gRe[(n - 1) * stride];
	GIm[(n - 1) * stride] = gIm[(n - 1) * stride];

	int j = 0; 	
	for(i = 0; i < n - 1; i++){
		GRe[i * stride] = gRe[j * stride]; 
		GIm[i * stride] = gIm[j * stride]; 
		int k = n/2; 
		while(k <= j){
			j -= k; 
			k /= 2; 
		} 
		j += k; 
	} 
	
	//Calculates FFT 

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
				double t1 = u1 * GRe[il * stride] - u2 * GIm[il * stride]; 
				double t2 = u1 * GIm[il * stride] - u2 * GRe[il * stride]; 
				GRe[il * stride] = GRe[i * stride] - t1; 
				GIm[il * stride] = GIm[i * stride] - t2; 
				GRe[i * stride] += t1; 
				GIm[i * stride] += t2; 
				

			} 
			double z = u1 * ca - u2 * sa; 
			u2 = u1 * sa + u2 * ca; 
			u1 = z; 

		} 
		sa = sqrt((1.0 - ca)/2.0); 
		if ( !inverse) sa =- sa; 
		ca = sqrt((1.0 + ca)/2.0); 
	} 
	if (!inverse){
		int i = 0; 
		for (i = 0; i < n; i++){
			GRe[i * stride] /= factor; 
			GIm[i * stride] /= factor; 
		} 
	} 



}  

void FFT2D(int w, int h, bool inverse, const double *gRe, const double *gIm, double *GRe, double *GIm){
	int val = w * h * 3;
	double Gr2[val]; 
	double Gi2[val]; 
	int y, c; 

	for (y = 0; y < h; y++){
		for ( c = 0; c < 3; c++){
			FFT(w, inverse, &gRe[y * w * 3 + c], &gIm[y * w * 3 + c], &Gr2[y * w * 3 + c],  &Gr2[y * w * 3 + c], 3, 1); 	
			printf("This works %d %d \n", y, c); 
		} 
	} 
	for (y = 0; y < h; y++){
		for ( c = 0; c < 3; c++){
			FFT(h, inverse, &Gr2[y * 3 + c], &Gi2[y * 3 + c], &GRe[y  * 3 + c],  &GIm[y * 3 + c], w * 3, inverse ? w : h); 	
		} 
	} 
} 
void load(double *fRe, double *fIm){
	unsigned char* data; 
	//unsigned int* mod; 

	data = stbi_load("test.png", &width, &height, &bpp, STBI_rgb_alpha); 

} 
void write(char *data, int width, int height, int comp, unsigned int *d){
//bpp = comp 
	int res; 
	res = stbi_write_png("output.png", width, height, comp, d, 0); 

} 

int main(){
	
	int i, j ; 

	for (i = 0; i < N; i++){
		for ( j = 0; j < N; j++){
			double x = 25 * sin(i /2.0);
			fRe[i][j][0] = x; 
			fIm[i][j][0] = x; 
			fIm[i][j][1] = x;  
			fRe[i][j][1] = x; 
			fRe[i][j][2] = x;  
			fIm[i][j][2] = x;  
		
		} 

	} 
	FFT2D(N, N, 0, fRe[0][0], fIm[0][0], FRe[0][0], FIm[0][0]);  
	load(fRe, fIm); 


	//int width, height, bpp; 
	//stbi_image_free(rgb); 
	return 0; 
} 
