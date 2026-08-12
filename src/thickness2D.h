/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#define MAX_ELEM_IN_BUCKET 200000
#define NUM_BUCKETS 400
#define MAXPATTERNS (16384*4)
#define MAXCLASSNUMBER MAXPATTERNS
#define MAXDIM 20
#ifndef UCHAR
#define UCHAR(c) ((unsigned char)(c))
#endif

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

struct element {
  int x;
  int y;
  int icur;
  int dcur;
  int idcur;
};

struct bucket {
  int num_elem;
  int *index_elem;
  int *index_l;
};

struct nodeDataNew {
  int id;
  int pclass;
  float d[MAXDIM];
  int row;
  int col;
  int slice;
};

int thickness2DYezzi(unsigned char* prototypes,int height, int width, float *maps, float** laplacefield,float** gradientx, float** gradienty, int num_it, float hx, float hy, unsigned char label_cortex,int debug);
int thickness2DYezzi_reverse(unsigned char* prototypes,int height, int width, float *maps, float** laplacefield,float** gradientx, float** gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug);
float** DToptimo(char* prototypes,int height, int width, int K, float** maps, int tipo_mapa);
float compute_mean_thickness2D(unsigned char *input, float *maps, int label_cortex, int height, int width, int *npoints, float *std);
