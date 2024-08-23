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
	double imag;
} Complex;

int main();

int img_process(RGB* img);

void print_prog(U32 cur_pos, U32 tgt);

int img_gain(RGB* img);

float fast_sqrt(float number);

RGB* load_bmp(const char* filename);

void save_bmp(const char* filename, RGB* img);

void dump_grey(U8* grey);

void fft1D(Complex* v, int n, Complex* tmp);

Complex* fft2D(U8* img, int width, int height);

int rgb2grey(RGB* img, U8* grey);


#endif
