/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include "utils.h"

#define MAX_ELEM_IN_BUCKET 200000
#define NUM_BUCKETS 400
#define MAXPATTERNS (16384 * 4)
#define MAXCLASSNUMBER MAXPATTERNS
#define MAXDIM 20
#ifndef UCHAR
#define UCHAR(c) ((unsigned char)(c))
#endif

#define round(x) ((x) >= 0 ? (long)((x) + 0.5) : (long)((x) - 0.5))

struct element {
  int x;
  int y;
  int icur;
  int dcur;
  int idcur;
};

struct bucket {
  int num_elem;
  int* index_elem;
  int* index_l;
};

struct nodeDataNew {
  int id;
  int pclass;
  float d[MAXDIM];
  int row;
  int col;
  int slice;
};

int thickness3DYezzi(unsigned char* prototypes, int height, int width, int depth, float* maps, float*** laplacefield, float*** gradientx, float*** gradienty, float*** gradientz, int num_it, float hx, float hy, float hz);
int thickness3DYezzi_reverse(unsigned char* prototypes, int height, int width, int depth, float* maps, float*** laplacefield, float*** gradientx, float*** gradienty, float*** gradientz, int num_it, float hx, float hy, float hz);
float** DToptimo(char* prototypes, int height, int width, int K, float** maps, int tipo_mapa);
float compute_mean_thickness(unsigned char* input, float* maps, int label_cortex, int boundary_l, int height, int width, int depth, float* sigma);
float compute_mean_thickness_volume(unsigned char* input, float* maps, int label_cortex, int height, int width, int depth, float* sigma);
