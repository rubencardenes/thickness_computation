/*(c) Ruben Cardenes Almeida, Valladolid, 8/11/2007 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "poisson2D.h"
#include "laplace2D.h"

/* Input values
255, o 3: Outside domain
1 Exterior boundary 
0 Interior boundary
2 Inside domain 
*/
int poisson2D(unsigned char* input, int height, int width, float** output, int iterations, float lambda, int reverse, float h) {
  int i, j, l;
  int sum = 0;
  init_laplace_field2D(input, height, width, output, reverse);

  /* Solve Poisson */
  for (l = 0; l < iterations; l++) {
    sum = 0;
    for (i = 0; i < height; i++) {
      for (j = 0; j < width; j++) {
        if (input[sum] == 2 && i != 0 && i != height - 1 && j != 0 && j != width - 1) {
          output[i][j] = 0.25 * (output[i - 1][j] + output[i + 1][j] + output[i][j - 1] + output[i][j + 1] - h * h);
          /*output[i][j] = output[i][j]+(lambda+1)*(0.25 *(output[i-1][j] + output[i+1][j] + output[i][j-1] + output[i][j+1]) - output[i][j]);*/
        }
        sum++;
      }
    }
  }

  return 0;
}


int minimos_locales2D(float** in, unsigned char* out, int height, int width) {
  int i, j, x, y, sum;
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      if (i == 0 || j == 0 || i == height - 1 || j == width - 1) continue;
      sum = 0;
      for (x = -1; x < 2; x++) {
        for (y = -1; y < 2; y++) {
          if (x == 0 && y == 0) continue;
          if (in[i + x][j + y] > in[i][j]) {
            sum++;
          }
        }
      }
      if (sum == 8) {
        out[i * width + j] = (int)in[i][j];
      } else {
        out[i * width + j] = 0;
      }
    }
  }

  return 0;
}
