// Factored discrete Fourier transform, or FFT, and its inverse iFFT
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
typedef struct cmplx{

double Re;
double Im;
	
}cplx;
cplx cexp(double n){

	cplx x;
	x.Re = cos(n);
	x.Im = sin(n);
	return x;
}
double cimag(cplx x){
	return x.Im;
}
double creal(cplx x){
	return x.Re;
}
void complexmult(cplx* x, cplx* y,cplx* result){

	result->Re = x->Re*y->Re - x->Im*y->Im;
	result->Im = x->Re*y->Im + x->Im*y->Re;
	
}
void _fft(cplx buf[], cplx out[], int n, int step)
{
	static double PI = 3.141592;


	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			cplx temp = cexp(-PI * i / n);
			cplx t;
			complexmult(&temp,&out[i+step],&t);
			buf[i / 2].Re = out[i].Re + t.Re;
			buf[i / 2].Im = out[i].Im + t.Im;
			buf[(i + n)/2].Re = out[i].Re - t.Re;
			buf[(i + n)/2].Im = out[i].Im - t.Im;
		}
	}
}
 
void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++){ 
		out[i].Re = buf[i].Re;
		out[i].Im = buf[i].Im;
	}
 
	_fft(buf, out, n, 1);
}
void show(const char * s, cplx buf[], int n) {
	printf("%s", s);
	for (int i = 0; i < n; i++)
		if (!cimag(buf[i]))
			printf("%f ", creal(buf[i]));
		else
			printf("(%f, %f) ", creal(buf[i]), cimag(buf[i]));
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
		reverse[n-i].Re = buf[i].Re;
		reverse[n-i].Im = buf[i].Im;
	}
	//printf("Reverse done\n");
	for (i = 0; i < n;i++)
	{
		//printf("Final i: %d\n",i);
		if (i < 1){
			(buf[i].Re /= n);
			(buf[i].Im /= n);
		}
		else
		{
			(buf[i].Re = reverse[i].Re/n);
			(buf[i].Im = reverse[i].Im/n);
		}
	}

	
}
void printcomp(cplx x,char *s){
	printf("%s Real: %f\t%s Imag: %f\n",s,x.Re,s,x.Im);}
void showf(char *s, cplx buf[], int n) {
	FILE *p;
	//printf("ShowF N: %d\n",n);
	p = fopen("/home/mert/Desktop/real-time-fft_c/result.txt","a");
	int i;
	fprintf(p,"%s: ", s);
	for (int i = 0; i < n; i++)
	{
			if (!cimag(buf[i]))
			fprintf(p,"%f ", creal(buf[i]));
		else
			fprintf(p,"(%f, %f) ", creal(buf[i]), cimag(buf[i]));
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
void randomize(cplx x[], int size, double range,double average,unsigned char imag){
	int i;
	for (i = 0; i < size;i++){
		x[i].Re = (rand() % (int)(2*range + 1)) - range + average;
		if (imag)
			x[i].Im = (rand() % (int)(2*range + 1)) - range + average;
		else
			x[i].Im = 0;
	}
	
}
void overlapsave(double *result,double *x,double *h,int L,int P, int M)
{
	int j;
	cplx fft_x[L];
	cplx ifft_x[L];
	cplx tmp_x[L];
	cplx fft_h[L];
	cplx x_zeropadded[P+L-(P%(L-M+1))];
	int i;
	int k;
	int maxincrement = P/(L- (M -1));
	if((P%(L-(M-1))) != 0)
		maxincrement++;
	for (i = 0; i < P+L-(P%(L-M+1)); i++)
	{
		if((i < (P + M - 1))&&(i >= (M-1))){
			x_zeropadded[i].Re= x[i-(M-1)];
			x_zeropadded[i].Im = 0;
		}
		else{
		 	x_zeropadded[i].Re = 0;
			x_zeropadded[i].Im = 0;
		}
	}
	for(i = 0; i < maxincrement ; i++)
	{
		for(j = 0; j < L; j++)
		{
			//printf("Second for: %d\n", (i*(L-(M-1))+j));
			tmp_x[j].Re = x_zeropadded[i*(L-(M-1))+j].Re;
			tmp_x[j].Im = x_zeropadded[i*(L-(M-1))+j].Im;

		}
		fft(tmp_x,L);
		for (k = 0; k< L; k++){
			fft_h[k].Re = (k < M) ? h[k] : 0.0;
			fft_h[k].Im = 0.0;
		}	
		fft(fft_h,L);
		for (k = 0;k < L; k++)
		{
			complexmult(&tmp_x[k],&fft_h[k],&ifft_x[k]);
		}
		ifft(ifft_x, L);
		for (k = M-1; k < L; k++)
			result[i*(L-(M-1))+(k - (M -1))] = creal(ifft_x[k]);

	}

}
int main()
{	
	int P = 1024; //Input size
	int M = 16; //filter length
	double *result;
	cplx x[P];
	cplx h[M];
	int seed = time(NULL);
    	srand(seed);
	int L = 128; // block size
	//x = malloc((P)*sizeof(cplx));
	//h = malloc((M)*sizeof(cplx));
	randomize(x,P,100,0,0); // Real signal decleration
	randomize(h,M,10,0,0);
	printf("Randomized\n");
	int i;
	overlapsave(result,x,h,L,P,M);
	printf("Overlapped and saved\n");
	char a[] = "Input";
	char b[] = "Filter";
	char c[] = "Result";
	showf(&a,x,P);
	showf(&b,h,M);
	showresult(&c,result,(P+M-1));
	return 0;
}
