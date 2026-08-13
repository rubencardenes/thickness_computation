/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#ifndef THICKNESS2D_H
#define THICKNESS2D_H

#include "utils.h"

int thickness2DYezzi(unsigned char* prototypes, int height, int width, float* maps, float** gradientx, float** gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug);
int thickness2DYezzi_reverse(unsigned char* prototypes, int height, int width, float* maps, float** gradientx, float** gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug);
int thickness2Dgradient(unsigned char* prototypes, int height, int width, float* input_maps, float* maps, float** gradientx, float** gradienty);
float compute_mean_thickness2D(unsigned char* input, float* maps, int label_cortex, int height, int width, int* npoints, float* std);

#endif /* THICKNESS2D_H */
