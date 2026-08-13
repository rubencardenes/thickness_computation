#ifndef LAPLACE2D_H
#define LAPLACE2D_H

#include "array_utils.h"
#include "utils.h"

/* Gradient definitions */
#define FDY(Data, iX, iY) (Data[(iX)][(iY) + 1] - Data[(iX)][(iY)])
#define BDY(Data, iX, iY) (Data[(iX)][(iY)] - Data[(iX)][(iY) - 1])

#define FDX(Data, iX, iY) (Data[(iX) + 1][(iY)] - Data[(iX)][(iY)])
#define BDX(Data, iX, iY) (Data[(iX)][(iY)] - Data[(iX) - 1][(iY)])

#define CDX(Data, iX, iY) ((Data[(iX) + 1][(iY)] - Data[(iX) - 1][(iY)]) * 0.5)
#define CDY(Data, iX, iY) ((Data[(iX)][(iY) + 1] - Data[(iX)][(iY) - 1]) * 0.5)


int new_compute_corners(unsigned short *input, int height, int width);
void init_laplace_field2D(const unsigned char *input, int height, int width, float **output, int reverse);
int laplace2D(unsigned char *input, int height, int width, float **output, int iterations, float lambda, int reverse);
int RelabelBoundary(unsigned char *domain, int height, int width);
int EdgeDetect(unsigned char *domain, int height, int width);
int normalize(float **gradientx, float **gradienty, int height, int width);
int iGradX(float **ppfData, float **ppfGradient, int numRowX, int numColY);
int iGradY(float **ppfData, float **ppfGradient, int numRowX, int numColY);
int maxcomponent2D(unsigned short *data, int height, int width, int label);
int sizefilter2D(unsigned short *data, int height, int width, int max_size, int oldlabel, int newlabel);
int floodfill(unsigned char *domain, int startindex, unsigned short oldlabel, unsigned char newlabel, int height, int width);
#endif /* LAPLACE2D_H */
