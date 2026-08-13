#ifndef LAPLACE3D_H
#define LAPLACE3D_H

#include "array_utils.h"
#include "utils.h"

/* 2D Gradient definitions */
#define FDY(Data, iX, iY) (Data[(iX)][(iY) + 1] - Data[(iX)][(iY)])
#define BDY(Data, iX, iY) (Data[(iX)][(iY)] - Data[(iX)][(iY) - 1])

#define FDX(Data, iX, iY) (Data[(iX) + 1][(iY)] - Data[(iX)][(iY)])
#define BDX(Data, iX, iY) (Data[(iX)][(iY)] - Data[(iX) - 1][(iY)])

#define CDX(Data, iX, iY) ((Data[(iX) + 1][(iY)] - Data[(iX) - 1][(iY)]) * 0.5)
#define CDY(Data, iX, iY) ((Data[(iX)][(iY) + 1] - Data[(iX)][(iY) - 1]) * 0.5)

/* 3D Gradient definitions */
#define FDY3(Data, iZ, iY, iX, hy) (Data[(iZ)][(iY) + 1][(iX)] - Data[(iZ)][(iY)][(iX)]) / hy
#define BDY3(Data, iZ, iY, iX, hy) (Data[(iZ)][(iY)][(iX)] - Data[(iZ)][(iY) - 1][(iX)]) / hy

#define FDZ3(Data, iZ, iY, iX, hz) (Data[(iZ) + 1][(iY)][(iX)] - Data[(iZ)][(iY)][(iX)]) / hz
#define BDZ3(Data, iZ, iY, iX, hz) (Data[(iZ)][(iY)][(iX)] - Data[(iZ) - 1][(iY)][(iX)]) / hz

#define FDX3(Data, iZ, iY, iX, hx) (Data[(iZ)][(iY)][(iX) + 1] - Data[(iZ)][(iY)][(iX)]) / hx
#define BDX3(Data, iZ, iY, iX, hx) (Data[(iZ)][(iY)][(iX)] - Data[(iZ)][(iY)][(iX) - 1]) / hx

#define CDZ3(Data, iZ, iY, iX, hz) ((Data[(iZ) + 1][(iY)][(iX)] - Data[(iZ) - 1][(iY)][(iX)]) / (hz * 2))
#define CDY3(Data, iZ, iY, iX, hy) ((Data[(iZ)][(iY) + 1][(iX)] - Data[(iZ)][(iY) - 1][(iX)]) / (hy * 2))
#define CDX3(Data, iZ, iY, iX, hx) ((Data[(iZ)][(iY)][(iX) + 1] - Data[(iZ)][(iY)][(iX) - 1]) / (hx * 2))

int laplace3D(unsigned char *input, int height, int width, int depth, float ***output, int iterations, float lambda);
int laplace3D_voxelsize(unsigned char *input, int height, int width, int depth, float ***output, int iterations, float hx, float hy, float hz, float lambda);
int RelabelBoundary3D(unsigned char *domain, int height, int width, int depth);
int EdgeDetect3D(unsigned char *domain, int height, int width, int depth);
int EdgeDetect3D_knee(unsigned char *domain, int height, int width, int depth);
int normalize3D(float ***gradientx, float ***gradienty, float ***gradientz, int height, int width, int depth);
int iGradX(float **ppfData, float **ppfGradient, int numRowX, int numColY);
int iGradY(float **ppfData, float **ppfGradient, int numRowX, int numColY);
int iGradX3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY, int numSlice, float hx);
int iGradY3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY, int numSlice, float hy);
int iGradZ3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY, int numSlice, float hz);
#endif /* LAPLACE3D_H */
