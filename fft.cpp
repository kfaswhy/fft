#include "fft.h"

BITMAPFILEHEADER fileHeader;
BITMAPINFOHEADER infoHeader;
int height = 0;
int width = 0;
int PaddingSize = 0;
BYTE* pad = NULL;

int iso = 1024;

int main()
{
	t0 = clock();

	RGB* img = NULL;

	char bmp_in[] = "C:/Work/Desktop/1.bmp";
	img = load_bmp(bmp_in);

	img_process(img);

	char bmp_proc[] = "C:/Work/Desktop/0_proc.bmp";
	save_bmp(bmp_proc, img);


	t1 = clock();
	U32 d_t = t1 - t0;
	LOG("sum time = %.3f.", (float)d_t / 1000);
	return 0;
}
 
int img_process(RGB* img)
{
	img_gain(img);

	U8* grey = (U8*)malloc(sizeof(U8) * height * width);

	rgb2grey(img, grey);

	dump_grey(grey);

	U8 img2[] = {
	1, 2, 3, 4,
	5, 6, 7, 8,
	9, 10, 11, 12,
	13, 14, 15, 16
	};
	Complex* data = fft2D(img2, 4, 4);

	return 0;
}

#if 1
void print_prog(U32 cur_pos, U32 tgt)
{
	end = clock();
	
	if ((end - start) >= 1000)
	{
		LOG("Processing: %d%%.", cur_pos * 100 / tgt);
		start = clock();
	}
}

float fast_sqrt(float number) {
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y = number;
	i = *(long*)&y;                       // 将 float 解释为 long 类型
	i = 0x5f3759df - (i >> 1);            // 魔术数字
	y = *(float*)&i;
	y = y * (threehalfs - (x2 * y * y));  // 近似值调整

	return 1.0f / y;
}

RGB* load_bmp(const char* filename)
{
	FILE* f_in = fopen(filename, "rb");

	fread(&fileHeader, sizeof(BITMAPFILEHEADER), 1, f_in);
	fread(&infoHeader, sizeof(BITMAPINFOHEADER), 1, f_in);

	height = infoHeader.biHeight;
	width = infoHeader.biWidth;
	int LineByteCnt = (((width * infoHeader.biBitCount) + 31) >> 5) << 2;
	//int ImageDataSize = LineByteCnt * height;
	PaddingSize = 4 - ((width * infoHeader.biBitCount) >> 3) & 3;
	pad = (BYTE*)malloc(sizeof(BYTE) * PaddingSize);
	RGB* img = (RGB*)malloc(sizeof(RGB) * height * width);

	if (infoHeader.biBitCount == 24) {
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int index = i * width + j;
				fread(&img[index], sizeof(RGB), 1, f_in);
			}
			if (PaddingSize != 0)
			{
				fread(pad, 1, PaddingSize, f_in);
			}
		}
	}
	else
	{
		printf("此程序不支持非24位图片");
		return NULL;
	}

	fclose(f_in);
	return img;
}

void save_bmp(const char* filename, RGB* img)
{
	FILE* f_out = fopen(filename, "wb");
	fwrite(&fileHeader, sizeof(fileHeader), 1, f_out);
	fwrite(&infoHeader, sizeof(infoHeader), 1, f_out);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++)
			fwrite(&img[i * width + j], sizeof(RGB), 1, f_out);
		if (PaddingSize != 0)
		{
			fwrite(pad, 1, PaddingSize, f_out);
		}
	}
	fclose(f_out);
	return;
}

void dump_grey(U8* grey)
{
	RGB* img = (RGB*)malloc(sizeof(RGB) * height * width);

	RGB* p_img = &img[0];
	U8* p_grey = &grey[0];

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			/* 像素运算 */
			p_img->r = *p_grey;
			p_img->g = *p_grey;
			p_img->b = *p_grey;
			p_img++;
			p_grey++;
			/* 像素运算结束 */
		}
	}

	char bmp_grey[] = "C:/Work/Desktop/grey.bmp";
	save_bmp(bmp_grey, img);
	free(img);
	return;
}
#endif


void fft1D(Complex* v, int n, Complex* tmp) {
	if (n > 1) {
		int k, m;
		Complex z, w, * vo, * ve;
		ve = tmp; vo = tmp + n / 2;
		for (k = 0; k < n / 2; k++) {
			ve[k] = v[2 * k];
			vo[k] = v[2 * k + 1];
		}
		fft1D(ve, n / 2, v);
		fft1D(vo, n / 2, v);
		for (m = 0; m < n / 2; m++) {
			w.real = cos(2 * PI * m / (double)n);
			w.imag = -sin(2 * PI * m / (double)n);
			z.real = w.real * vo[m].real - w.imag * vo[m].imag;
			z.imag = w.real * vo[m].imag + w.imag * vo[m].real;
			v[m].real = ve[m].real + z.real;
			v[m].imag = ve[m].imag + z.imag;
			v[m + n / 2].real = ve[m].real - z.real;
			v[m + n / 2].imag = ve[m].imag - z.imag;
		}
	}
}

Complex* fft2D(U8* img, int width, int height) {
	int i, j;
	Complex* data = (Complex*)malloc(width * height * sizeof(Complex));
	Complex* tmp = (Complex*)malloc(width * sizeof(Complex));

	// Convert image data to complex numbers
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			data[i * width + j].real = img[i * width + j];
			data[i * width + j].imag = 0.0;
		}
	}

	// Perform FFT on rows
	for (i = 0; i < height; i++) {
		fft1D(&data[i * width], width, tmp);
	}

	// Perform FFT on columns
	tmp = (Complex*)realloc(tmp, height * sizeof(Complex));
	for (j = 0; j < width; j++) {
		for (i = 0; i < height; i++) {
			tmp[i] = data[i * width + j];
		}
		fft1D(tmp, height, data);
		for (i = 0; i < height; i++) {
			data[i * width + j] = tmp[i];
		}
	}

	// Print the result (for demonstration purposes)
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			printf("(%.02f, %.02f) ", data[i * width + j].real, data[i * width + j].imag);
		}
		printf("\n");
	}

	free(tmp);
	return data;
}



int rgb2grey(RGB* img, U8* grey)
{
	RGB* p_img = &img[0];
	U8* p_grey = &grey[0];

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			/* 像素运算 */
			*p_grey = (p_img->g + p_img->r + p_img->b)/3;
			p_img++;
			p_grey++;
			/* 像素运算结束 */
		}
	}
	return 0;
}


int img_gain(RGB* img)
{
	if (iso == 1024)
	{
		return 0;
	}
	
	RGB* p_img = &img[0];
	
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			/* 像素运算 */
			p_img->r = clp_range(0, (p_img->r * iso + 512) >> 10, U8MAX);
			p_img->g = clp_range(0, (p_img->g * iso + 512) >> 10, U8MAX);
			p_img->b = clp_range(0, (p_img->b * iso + 512) >> 10, U8MAX);

			p_img++;
			/* 像素运算结束 */
		}
	}
	return 0;
}