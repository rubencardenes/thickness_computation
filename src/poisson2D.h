#ifndef POISSON2D_H
#define POISSON2D_H


int poisson2D(unsigned char* input, int height, int width, float** output, int iterations, float lambda, int reverse, float h, int cortex_label);
#endif /* POISSON2D_H */
