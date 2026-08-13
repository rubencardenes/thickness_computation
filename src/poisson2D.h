#ifndef POISSON2D_H
#define POISSON2D_H


int poisson2D(unsigned char* input, int height, int width, float** output, int iterations, float lambda, int reverse, float h);
/* Marks the strict local minima of `in` (8-neighbourhood, excluding the
   outermost 1-pixel border) into `out`, storing the truncated field value
   there and 0 everywhere else. */
int minimos_locales2D(float** in, unsigned char* out, int height, int width);
#endif /* POISSON2D_H */
