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

	//char bmp_proc[] = "C:/Work/Desktop/0_proc.bmp";
	//save_bmp(bmp_proc, img);


	t1 = clock();
	U32 d_t = t1 - t0;
	LOG("sum time = %.3f.", (float)d_t / 1000);
	return 0;
}
 
int img_process(RGB* img)
{
	img_gain(img);

	double* grey = (double*)malloc(sizeof(double) * height * width);

	rgb2grey(img, grey);

	dump_grey(grey); 


	Complex* fd = (Complex*)malloc(sizeof(Complex) * height * width);
	//FFT2D(data, dst, 4, 4);
	FFT2D(grey, fd, width, height);
	RGB* mag = (RGB*)malloc(sizeof(RGB) * height * width);
	RGB* deg = (RGB*)malloc(sizeof(RGB) * height * width);


	//filter_complex(fd);
	
	dump_complex(fd, mag, deg, 1);
	char bmp_mag[] = "C:/Work/Desktop/2 mag.bmp";
	save_bmp(bmp_mag, mag);
	char bmp_deg[] = "C:/Work/Desktop/3 deg.bmp";
	save_bmp(bmp_deg, deg);


	Complex* ifd = (Complex*)malloc(sizeof(Complex) * height * width);
	IFFT2D(fd, ifd, width, height);
	//Complex* data = fft2D(img2, 4, 4);


	RGB* proc = (RGB*)malloc(sizeof(RGB) * height * width);
	RGB* proc2 = (RGB*)malloc(sizeof(RGB) * height * width);
	dump_complex(ifd, proc, proc2, 0);
	char bmp_prc[] = "C:/Work/Desktop/4 proc.bmp";
	save_bmp(bmp_prc, proc);


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

void dump_grey(double* grey)
{
	RGB* img = (RGB*)malloc(sizeof(RGB) * height * width);

	RGB* p_img = &img[0];
	double* p_grey = &grey[0];

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			/* 像素运算 */
			*p_grey = clp_range(0, *p_grey, 255);
			p_img->r = (U8)*p_grey;
			p_img->g = (U8)*p_grey;
			p_img->b = (U8)*p_grey;
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

void dump_complex(Complex* data, RGB* mag, RGB* deg, U8 is_revert)
{
	RGB* mag_revert = (RGB*)malloc(sizeof(RGB) * height * width);

	RGB* p_mag = &mag[0];
	RGB* p_deg = &deg[0];
	Complex* p_data = &data[0];

	double mag2_max = 0.0, mag2_min = 0.0, mag2 = 0.0;
	double degree = 0.0;

	//计算极值
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			mag2 = p_data->real * p_data->real + p_data->imagin * p_data->imagin;
			mag2_max = calc_max(mag2, mag2_max);

			//LOG("%03f, %03f.", p_data->real, p_data->imagin);
			p_data++;
		}
	}

	LOG("mag2max = %3f.", mag2_max);

	p_data = &data[0];
	//double range = 255.0 / (mag2_max - mag2_min);
	double tmp = 0;


	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			mag2 = p_data->real * p_data->real + p_data->imagin * p_data->imagin;
			mag2 = mag2 * 65025 / mag2_max;

			p_mag->r = (U8)clp_range(0, fast_sqrt(mag2), 255);
			p_mag->g = p_mag->r;
			p_mag->b = p_mag->r;
			//tmp = p_data->mag - mag_min;
			//tmp *= range;
			//LOG("%f.", tmp);

			if (p_data->real == 0)
			{
				degree = 0;
			}
			else
			{
				degree = atan(p_data->imagin / p_data->real) * 180 / PI;
				while (degree < 0)
				{
					degree += 360;
				}
				while (degree >= 360)
				{
					degree -= 360;
				}
			}

			p_deg->r = (U8)clp_range(0, degree * 255 / 360, 255);
			p_deg->g = p_deg->r;
			p_deg->b = p_deg->r;

			p_data++;
			p_mag++;
			p_deg++;
		}
	}

	if (is_revert)
	{
		for (int i = 0; i < height / 2; i++)
		{
			for (int j = 0; j < width / 2; j++)
			{
				int index_mag = i * width + j;
				int index_revert = (i + height / 2) * width + (j + width / 2);
				mag_revert[index_revert] = mag[index_mag];
			}
		}

		for (int i = 0; i < height / 2; i++)
		{
			for (int j = width / 2; j < width; j++)
			{
				int index_mag = i * width + j;
				int index_revert = (i + height / 2) * width + (j - width / 2);
				mag_revert[index_revert] = mag[index_mag];
			}
		}

		for (int i = height / 2; i < height; i++)
		{
			for (int j = 0; j < width / 2; j++)
			{
				int index_mag = i * width + j;
				int index_revert = (i - height / 2) * width + (j + width / 2);
				mag_revert[index_revert] = mag[index_mag];
			}
		}

		for (int i = height / 2; i < height; i++)
		{
			for (int j = width / 2; j < width; j++)
			{
				int index_mag = i * width + j;
				int index_revert = (i - height / 2) * width + (j - width / 2);
				mag_revert[index_revert] = mag[index_mag];
			}
		}


		memcpy(mag, mag_revert, sizeof(RGB) * width * height);
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int index = i * width + j;
				mag_revert[index] = mag[index];
				mag[index].r = clp_range(0, log((double)mag_revert[index].r + 1) * 255 / 5.55, 255);
				mag[index].g = mag[index].r;
				mag[index].b = mag[index].r;
			}
		}
	}
	free(mag_revert);
	return;
}
#endif


void filter_complex(Complex* img)
{
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int tmpi = i < (height / 2) ? i : (height - i);
			int tmpj = j < (width / 2) ? j : (width - j);
			int index = i * width + j;
			int tmp = tmpi * tmpi + tmpj * tmpj;
			if (tmp > (width * width / 4 - 1000)&&
				tmp <  (width * width / 4 + 1000))
			{
				
				img[index].real = 255;
				img[index].imagin = 0;
			}
			else if (tmp > (width * width / 16 - 1000) &&
				tmp < (width * width / 16 + 1000))
			{

				img[index].real = 150;
				img[index].imagin = 0;
			}
			else
			{
				img[index].real = 0;
				img[index].imagin = 0;
			}

			

		}
	}

	//for (int i = 44; i < 84; i++) {
	//	for (int j = 408; j < 488; j++) {
	//		int index = i * width + j;
	//		img[index].real = 0;
	//		img[index].imagin = 0;
	//	}
	//}
}




int isBase2(int size_n) {
	int k = size_n;
	int z = 0;
	while (k /= 2) {
		z++;
	}
	k = z;
	if (size_n != (1 << k))
		return -1;
	else
		return k;
}
////////////////////////////////////////////////////////////////////
//复数基本运算
///////////////////////////////////////////////////////////////////
void Add_Complex(Complex* src1, Complex* src2, Complex* dst) {
	dst->imagin = src1->imagin + src2->imagin;
	dst->real = src1->real + src2->real;
}
void Sub_Complex(Complex* src1, Complex* src2, Complex* dst) {
	dst->imagin = src1->imagin - src2->imagin;
	dst->real = src1->real - src2->real;
}
void Multy_Complex(Complex* src1, Complex* src2, Complex* dst) {
	double r1 = 0.0, r2 = 0.0;
	double i1 = 0.0, i2 = 0.0;
	r1 = src1->real;
	r2 = src2->real;
	i1 = src1->imagin;
	i2 = src2->imagin;
	dst->imagin = r1 * i2 + r2 * i1;
	dst->real = r1 * r2 - i1 * i2;
}
void Copy_Complex(Complex* src, Complex* dst) {
	dst->real = src->real;
	dst->imagin = src->imagin;
}
void Show_Complex(Complex* src, int size_n) {
	if (size_n == 1) {
		if (src->imagin >= 0.0)
			printf("%lf+%lfj  ", src->real, src->imagin);
		else
			printf("%lf%lfj  ", src->real, src->imagin);

	}
	else if (size_n > 1) {
		for (int i = 0; i < size_n; i++)
			if (src[i].imagin >= 0.0) {
				printf("%lf+%lfj  ", src[i].real, src[i].imagin);
			}
			else
				printf("%lf%lfj  ", src[i].real, src[i].imagin);



	}


}
////////////////////////////////////////////////////////////////////
//计算WN
///////////////////////////////////////////////////////////////////
void getWN(double n, double size_n, Complex* dst) {
	double x = 2.0 * M_PI * n / size_n;
	dst->imagin = -sin(x);
	dst->real = cos(x);
}
////////////////////////////////////////////////////////////////////
//随机初始化输入
///////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
//标准DFT
///////////////////////////////////////////////////////////////////
void DFT(double* src, Complex* dst, int size) {

	for (int m = 0; m < size; m++) {
		double real = 0.0;
		double imagin = 0.0;
		for (int n = 0; n < size; n++) {
			double x = M_PI * 2 * m * n;
			real += src[n] * cos(x / size);
			imagin += src[n] * (-sin(x / size));

		}
		dst[m].imagin = imagin;
		dst[m].real = real;

	}
}
////////////////////////////////////////////////////////////////////
//IDT，复原傅里叶
///////////////////////////////////////////////////////////////////
void IDFT(Complex* src, Complex* dst, int size) {
	for (int m = 0; m < size; m++) {
		double real = 0.0;
		double imagin = 0.0;
		for (int n = 0; n < size; n++) {
			double x = M_PI * 2 * m * n / size;
			real += src[n].real * cos(x) - src[n].imagin * sin(x);
			imagin += src[n].real * sin(x) + src[n].imagin * cos(x);

		}
		real /= size;
		imagin /= size;
		if (dst != NULL) {
			dst[m].real = real;
			dst[m].imagin = imagin;
		}
	}


}
////////////////////////////////////////////////////////////////////
//FFT前，对输入数据重新排序
///////////////////////////////////////////////////////////////////
int FFTReal_remap(double* src, int size_n) {

	if (size_n == 1)
		return 0;
	double* temp = (double*)malloc(sizeof(double) * size_n);
	for (int i = 0; i < size_n; i++)
		if (i % 2 == 0)
			temp[i / 2] = src[i];
		else
			temp[(size_n + i) / 2] = src[i];
	for (int i = 0; i < size_n; i++)
		src[i] = temp[i];
	free(temp);
	FFTReal_remap(src, size_n / 2);
	FFTReal_remap(src + size_n / 2, size_n / 2);
	return 1;


}
int FFTComplex_remap(Complex* src, int size_n) {

	if (size_n == 1)
		return 0;
	Complex* temp = (Complex*)malloc(sizeof(Complex) * size_n);
	for (int i = 0; i < size_n; i++)
		if (i % 2 == 0)
			Copy_Complex(&src[i], &(temp[i / 2]));
		else
			Copy_Complex(&(src[i]), &(temp[(size_n + i) / 2]));
	for (int i = 0; i < size_n; i++)
		Copy_Complex(&(temp[i]), &(src[i]));
	free(temp);
	FFTComplex_remap(src, size_n / 2);
	FFTComplex_remap(src + size_n / 2, size_n / 2);
	return 1;


}


int DFT2D(double* src, Complex* dst, int size_w, int size_h) {
	for (int u = 0; u < size_w; u++) {
		for (int v = 0; v < size_h; v++) {
			double real = 0.0;
			double imagin = 0.0;
			for (int i = 0; i < size_w; i++) {
				for (int j = 0; j < size_h; j++) {
					double I = src[i * size_w + j];
					double x = M_PI * 2 * ((double)i * u / (double)size_w + (double)j * v / (double)size_h);
					real += cos(x) * I;
					imagin += -sin(x) * I;

				}
			}
			dst[u * size_w + v].real = real;
			dst[u * size_w + v].imagin = imagin;

		}

	}
	return 0;
}
/*

 */
int IDFT2D(Complex* src, Complex* dst, int size_w, int size_h) {
	for (int i = 0; i < size_w; i++) {
		for (int j = 0; j < size_h; j++) {
			double real = 0.0;
			double imagin = 0.0;
			for (int u = 0; u < size_w; u++) {
				for (int v = 0; v < size_h; v++) {
					double R = src[u * size_w + v].real;
					double I = src[u * size_w + v].imagin;
					double x = M_PI * 2 * ((double)i * u / (double)size_w + (double)j * v / (double)size_h);
					real += R * cos(x) - I * sin(x);
					imagin += I * cos(x) + R * sin(x);

				}
			}
			dst[i * size_w + j].real = (1. / (size_w * size_h)) * real;
			dst[i * size_w + j].imagin = (1. / (size_w * size_h)) * imagin;

		}
	}
	return 0;
}
/*



 */
void ColumnVector(Complex* src, Complex* dst, int size_w, int size_h) {
	for (int i = 0; i < size_h; i++)
		Copy_Complex(&src[size_w * i], &dst[i]);

}
/*

 */
void IColumnVector(Complex* src, Complex* dst, int size_w, int size_h) {
	for (int i = 0; i < size_h; i++)
		Copy_Complex(&src[i], &dst[size_w * i]);

}
/*
 */
int FFT2D(double* src, Complex* dst, int size_w, int size_h) {
	if (isBase2(size_w) == -1 || isBase2(size_h) == -1)
		exit(0);
	Complex* temp = (Complex*)malloc(sizeof(Complex) * size_h * size_w);
	if (temp == NULL)
		return -1;
	for (int i = 0; i < size_h; i++) {
		RealFFT(&src[size_w * i], &temp[size_w * i], size_w);
	}

	Complex* Column = (Complex*)malloc(sizeof(Complex) * size_h);
	if (Column == NULL)
		return -1;
	for (int i = 0; i < size_w; i++) {
		ColumnVector(&temp[i], Column, size_w, size_h);
		FFT(Column, Column, size_h);
		IColumnVector(Column, &temp[i], size_w, size_h);

	}



	for (int i = 0; i < size_h * size_w; i++)
		Copy_Complex(&temp[i], &dst[i]);
	free(temp);
	free(Column);
	return 0;
}

////////////////////////////////////////////////////////////////////
//FFT公式
///////////////////////////////////////////////////////////////////
void FFT(Complex* src, Complex* dst, int size_n) {

	int k = size_n;
	int z = 0;
	while (k /= 2) {
		z++;
	}
	k = z;
	if (size_n != (1 << k))
		exit(0);
	Complex* src_com = (Complex*)malloc(sizeof(Complex) * size_n);
	if (src_com == NULL)
		exit(0);
	for (int i = 0; i < size_n; i++) {
		Copy_Complex(&src[i], &src_com[i]);
	}
	FFTComplex_remap(src_com, size_n);
	for (int i = 0; i < k; i++) {
		z = 0;
		for (int j = 0; j < size_n; j++) {
			if ((j / (1 << i)) % 2 == 1) {
				Complex wn;
				getWN(z, size_n, &wn);
				Multy_Complex(&src_com[j], &wn, &src_com[j]);
				z += 1 << (k - i - 1);
				Complex temp;
				int neighbour = j - (1 << (i));
				temp.real = src_com[neighbour].real;
				temp.imagin = src_com[neighbour].imagin;
				Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
				Sub_Complex(&temp, &src_com[j], &src_com[j]);
			}
			else
				z = 0;
		}

	}


	for (int i = 0; i < size_n; i++) {
		Copy_Complex(&src_com[i], &dst[i]);
	}
	free(src_com);


}
void RealFFT(double* src, Complex* dst, int size_n) {


	int k = size_n;
	int z = 0;
	while (k /= 2) {
		z++;
	}
	k = z;
	if (size_n != (1 << k))
		exit(0);
	Complex* src_com = (Complex*)malloc(sizeof(Complex) * size_n);
	if (src_com == NULL)
		exit(0);
	for (int i = 0; i < size_n; i++) {
		src_com[i].real = src[i];
		src_com[i].imagin = 0;
	}
	FFTComplex_remap(src_com, size_n);
	for (int i = 0; i < k; i++) {
		z = 0;
		for (int j = 0; j < size_n; j++) {
			if ((j / (1 << i)) % 2 == 1) {
				Complex wn;
				getWN(z, size_n, &wn);
				Multy_Complex(&src_com[j], &wn, &src_com[j]);
				z += 1 << (k - i - 1);
				Complex temp;
				int neighbour = j - (1 << (i));
				temp.real = src_com[neighbour].real;
				temp.imagin = src_com[neighbour].imagin;
				Add_Complex(&temp, &src_com[j], &src_com[neighbour]);
				Sub_Complex(&temp, &src_com[j], &src_com[j]);
			}
			else
				z = 0;
		}

	}


	for (int i = 0; i < size_n; i++) {
		Copy_Complex(&src_com[i], &dst[i]);
	}
	free(src_com);

}

void IFFT(Complex* src, Complex* dst, int size_n) {
	for (int i = 0; i < size_n; i++)
		src[i].imagin = -src[i].imagin;
	FFTComplex_remap(src, size_n);
	int z, k;
	if ((z = isBase2(size_n)) != -1)
		k = isBase2(size_n);
	else
		exit(0);
	for (int i = 0; i < k; i++) {
		z = 0;
		for (int j = 0; j < size_n; j++) {
			if ((j / (1 << i)) % 2 == 1) {
				Complex wn;
				getWN(z, size_n, &wn);
				Multy_Complex(&src[j], &wn, &src[j]);
				z += 1 << (k - i - 1);
				Complex temp;
				int neighbour = j - (1 << (i));
				Copy_Complex(&src[neighbour], &temp);
				Add_Complex(&temp, &src[j], &src[neighbour]);
				Sub_Complex(&temp, &src[j], &src[j]);
			}
			else
				z = 0;
		}

	}
	for (int i = 0; i < size_n; i++) {
		dst[i].imagin = (1. / size_n) * src[i].imagin;
		dst[i].real = (1. / size_n) * src[i].real;
	}
}

/*
 */
int IFFT2D(Complex* src, Complex* dst, int size_w, int size_h) {

	if (isBase2(size_w) == -1 || isBase2(size_h) == -1)
		exit(0);

	Complex* temp = (Complex*)malloc(sizeof(Complex) * size_h * size_w);
	if (temp == NULL)
		return -1;
	Complex* Column = (Complex*)malloc(sizeof(Complex) * size_h);
	if (Column == NULL)
		return -1;

	for (int i = 0; i < size_w; i++) {
		ColumnVector(&src[i], Column, size_w, size_h);
		IFFT(Column, Column, size_h);
		IColumnVector(Column, &src[i], size_w, size_h);

	}
	for (int i = 0; i < size_w * size_h; i++)
		src[i].imagin = -src[i].imagin;
	for (int i = 0; i < size_h; i++) {
		IFFT(&src[size_w * i], &temp[size_w * i], size_w);

	}


	for (int i = 0; i < size_h * size_w; i++)
	{
		Copy_Complex(&temp[i], &dst[i]);
		//printf("%.3f, %.3f\n", dst[i].real, dst[i].imagin);
	}

	free(temp);
	free(Column);
	return 0;

}



int rgb2grey(RGB* img, double* grey)
{
	RGB* p_img = &img[0];
	double* p_grey = &grey[0];

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