#ifndef _DEFOG_H_
#define _DEFOG_H_

#include <iostream>
#include <windows.h> 
#include <time.h>
#include <math.h>

#define U64 unsigned long long
#define U32 unsigned int
#define S32 int
#define U16 unsigned short
#define U8 unsigned char 

constexpr auto U16MAX = (0xFFFF);
constexpr auto U8MAX = (255);
constexpr auto U8MIN = (0);
constexpr auto PI = 3.1415926;
constexpr auto M_PI = 3.1415926;

#define calc_min(a,b) ((a)>(b)?(b):(a))
#define calc_max(a,b) ((a)<(b)?(b):(a))
#define calc_abs(a) ((a)>0?(a):(-a))
#define clp_range(min,x,max) calc_min(calc_max((x), (min)), (max))

#define LOG(...) printf("%s [%d]: ", __FUNCTION__, __LINE__);printf(__VA_ARGS__);printf("\n");


U32 start = clock();
U32 end;
U8 prog_print = 1;

U32 t0;
U32 t1;

typedef struct _RGB
{
	BYTE b;
	BYTE g;
	BYTE r;
}RGB;


typedef struct {
	double real;
	double imagin;
} Complex;

int main();

int img_process(RGB* img);

void print_prog(U32 cur_pos, U32 tgt);

int img_gain(RGB* img);

float fast_sqrt(float number);

RGB* load_bmp(const char* filename);

void save_bmp(const char* filename, RGB* img);

void dump_grey(double* grey);

void dump_complex(Complex* data, RGB* mag, RGB* deg, U8 is_revert);
void filter_complex(Complex* img);
int isBase2(int size_n);

void Add_Complex(Complex* src1, Complex* src2, Complex* dst);

void Sub_Complex(Complex* src1, Complex* src2, Complex* dst);

void Multy_Complex(Complex* src1, Complex* src2, Complex* dst);

void Copy_Complex(Complex* src, Complex* dst);

void Show_Complex(Complex* src, int size_n);

void getWN(double n, double size_n, Complex* dst);

void DFT(double* src, Complex* dst, int size);

void IDFT(Complex* src, Complex* dst, int size);

int FFTReal_remap(double* src, int size_n);

int FFTComplex_remap(Complex* src, int size_n);

int DFT2D(double* src, Complex* dst, int size_w, int size_h);

int IDFT2D(Complex* src, Complex* dst, int size_w, int size_h);

void ColumnVector(Complex* src, Complex* dst, int size_w, int size_h);

void IColumnVector(Complex* src, Complex* dst, int size_w, int size_h);

int FFT2D(double* src, Complex* dst, int size_w, int size_h);

void FFT(Complex* src, Complex* dst, int size_n);

void RealFFT(double* src, Complex* dst, int size_n);

void IFFT(Complex* src, Complex* dst, int size_n);

int IFFT2D(Complex* src, Complex* dst, int size_w, int size_h);


int rgb2grey(RGB* img, double* grey);


#endif
