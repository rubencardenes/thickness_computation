/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#ifndef THICKNESS3D_H
#define THICKNESS3D_H

#include "utils.h"

int thickness3DYezzi(unsigned char* prototypes, int height, int width, int depth, float* maps, float*** gradientx, float*** gradienty, float*** gradientz, int num_it, float hx, float hy, float hz);
int thickness3DYezzi_reverse(unsigned char* prototypes, int height, int width, int depth, float* maps, float*** gradientx, float*** gradienty, float*** gradientz, int num_it, float hx, float hy, float hz);
float compute_mean_thickness(unsigned char* input, float* maps, int label_cortex, int boundary_l, int height, int width, int depth, float* sigma);
float compute_mean_thickness_volume(unsigned char* input, float* maps, int label_cortex, int height, int width, int depth, float* sigma);

#endif /* THICKNESS3D_H */
