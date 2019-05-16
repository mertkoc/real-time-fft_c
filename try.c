// Factored discrete Fourier transform, or FFT, and its inverse iFFT
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#define PI 3.141592
typedef double complex cplx;
cplx *x;
cplx *h;
void _fft(cplx buf[], cplx out[], int n, int step)
{


	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
 
void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++) out[i] = buf[i];
 
	_fft(buf, out, n, 1);
}
void show(const char * s, cplx buf[], int n) {
	printf("%s", s);
	for (int i = 0; i < n; i++)
		if (!cimag(buf[i]))
			printf("%g ", creal(buf[i]));
		else
			printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}
void ifft(cplx buf[], int n)
{
	int i;
	//printf("N = %d\n",n);
	cplx reverse[n];
	fft(buf,n);
	//printf("FFT done\n");
	for (i = 1; i < n;i++)
	{
		reverse[n-i] = buf[i];
	}
	//printf("Reverse done\n");
	for (i = 0; i < n;i++)
	{
		//printf("Final i: %d\n",i);
		i < 1 ? (buf[i] /= n) : (buf[i] = reverse[i]/n);
	}

	
}
void showf(char *s, cplx buf[], int n) {
	FILE *p;
	//printf("ShowF N: %d\n",n);
	p = fopen("/home/mert/Desktop/real-time-fft_c/result.txt","a");
	int i;
	fprintf(p,"%s: ", s);
	for (int i = 0; i < n; i++)
	{
			if (!cimag(buf[i]))
			fprintf(p,"%g ", creal(buf[i]));
		else
			fprintf(p,"(%g, %g) ", creal(buf[i]), cimag(buf[i]));
	}
	fprintf(p,"\n");
	fclose(p);
}
void showresult(char *s, double *result,int n)
{
	FILE *p;
	p = fopen("/home/mert/Desktop/real-time-fft_c/result.txt","a");
	int i;
	fprintf(p,"%s: ", s);
	for (int i = 0; i < n; i++)
	{
		fprintf(p,"%f\t",result[i]);
	}
	fprintf(p,"\n");
	fclose(p);
}
double *overlapsave(cplx x[],cplx h[],int L,int P, int M)
{
	int j;
	cplx fft_x[L];
	cplx ifft_x[L];
	cplx tmp_x[L];
	cplx fft_h[L];
		
	
	cplx x_zeropadded[P+M];
	int i;
	int k;
	int maxincrement = P/(L- (M -1));
	double result [P+M-1];
    for(i=0;i<P;i++)
    {
		//printf("%d ",i);
		x[i] = rand() % 201 -100;
	}
	for (i = 0; i<M; i++)
	{
		h[i] = rand() % 10;
		
	}
	//printf("Filter yaratti\n");
	for (i = 0; i < P+M-1;i++)
	{
		if( i < M -1)
			x_zeropadded[i] = 0;
		else
			x_zeropadded[i] = x[i-(M-1)];
	}
	//printf("Zeropadded\n");
	for(i = 0; i < maxincrement ; i++)
	{
		//printf("First for: %d\n", i);
		for(j = 0; j < L; j++)
		{
			//printf("Second for: %d\n", (i*(L-(M-1))+j));
			tmp_x[j] = x_zeropadded[i*(L-(M-1))+j];
			//printf("Cikti\n");
			
		}
		fft(tmp_x,L);
		for (k = 0; k< L; k++)
			fft_h[k] = h[k];
		fft(fft_h,L);
		for (k = 0;k < L; k++)
			ifft_x[k] = tmp_x[k]*fft_h[k];
		//printf("IFFT\n");
		ifft(ifft_x, L);
		for (k = M-1; k < L; k++)
			result[i*(L-(M-1))+(k - (M -1))] = creal(ifft_x[k]);
		
	}
	for(i = 0; i < (P+M-1) ; i++)
		printf("%f\t",result[i]);
	printf("\n");
	return result;

}
int main()
{

	int P = 1024; //Input size
	int M = 16; //filter length
	double *result;
	int seed = time(NULL);
    srand(seed);
	int L = 128; // block size
	result = malloc((P+M-1)*sizeof(cplx));
	x = malloc((P)*sizeof(cplx));
	h = malloc((M)*sizeof(cplx));
	//printf("Sizeof x = %ld\n",(sizeof(x)/sizeof(x[0])));
	//printf("Overlapsave basladi\n");
	result = overlapsave(x,h,L,P,M);
	
	char a[] = "Input";
	char b[] = "Filter";
	//char c[] = "Result";
	showf(&a,x,P);
	showf(&b,h,M);
	//showresult(&c,result,(P+M-1));
	
	
	return 0;
}
