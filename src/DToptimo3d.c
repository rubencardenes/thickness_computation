/* Copyright (c) Ruben Cardenes Almeida 22/03/2002 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "utils.h"

#define MAX_ELEM_IN_BUCKET 1000000
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
  int z;
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

static int numelembucket[NUM_BUCKETS];

float distance3d_voxelsize(int x1, int y1, int z1, int x2, int y2, int z2, float hx, float hy, float hz) {
  return sqrt((x1 - x2) * (x1 - x2) * hx * hx + (y1 - y2) * (y1 - y2) * hy * hy + (z1 - z2) * (z1 - z2) * hz * hz);
}

int assign(struct element *Element_p, int l, int mapindex, int d, float **maps, int K) {
  maps[(*Element_p).icur][mapindex] = (char)l;
  /* si queremos mapas de distancias */
  maps[K][mapindex] = d;

  (*Element_p).icur++;
  if ((*Element_p).dcur != d) {
    (*Element_p).idcur = (*Element_p).icur;
  }
  (*Element_p).dcur = d;
  return 0;
}

/* Pushes the up-to-6 face neighbours of Element_p that are still domain, and
   whose distance to prototype `l` lands in the current or next bucket, onto
   those buckets. Also keeps maps[2] as the running minimum real distance.

   This was three near-identical 38-line blocks, one per axis; the only thing
   that varied was which of x/y/z carried the -1/+1 step. `offset` lists the
   neighbours in the order those blocks visited them, which the bucket contents
   depend on. */
int propagate3d(struct element Element_p, int l, int height, int width, int depth, struct element *proto, struct bucket *Bucket, float **maps, int dcur, struct element *Element, int K, int domain_label, char *prototypes, float hx, float hy, float hz) {
  static const int offset[6][3] = {
      {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
  int new_mapindex, dist, siguiente, n;
  int x, y, z, nx, ny, nz;
  int icuraux;
  float dreal;

  for (n = 0; n < 6; n++) {
    x = offset[n][0];
    y = offset[n][1];
    z = offset[n][2];
    nx = Element_p.x + x;
    ny = Element_p.y + y;
    nz = Element_p.z + z;

    if (nx < 0 || ny < 0 || nz < 0 || nx >= height || ny >= width || nz >= depth) continue;

    /* Only step away from the prototype, never back towards it. */
    if (!(abs(Element_p.x - proto[l].x) < abs(nx - proto[l].x) ||
          abs(Element_p.y - proto[l].y) < abs(ny - proto[l].y) ||
          abs(Element_p.z - proto[l].z) < abs(nz - proto[l].z))) continue;

    new_mapindex = mapIndex3D(nx, ny, nz, height, width, depth);
    if (prototypes[new_mapindex] != domain_label) continue;

    icuraux = Element[new_mapindex].icur > 0 ? Element[new_mapindex].icur - 1 : 0;
    if (!(Element[new_mapindex].icur < K && maps[icuraux][new_mapindex] != l)) continue;

    dreal = distance3d_voxelsize(nx, ny, nz, proto[l].x, proto[l].y, proto[l].z, hx, hy, hz);
    dist = round(dreal);
    if (dist > dcur + 1 || dist < dcur) continue;

    /* put p,l in Bucket d=dist */
    siguiente = Bucket[dist].num_elem;
    if (siguiente >= MAX_ELEM_IN_BUCKET) {
      printf("Excedidos num elem in bucket %d\n", MAX_ELEM_IN_BUCKET);
      return 1;
    }
    Bucket[dist].index_elem[siguiente] = new_mapindex;
    Bucket[dist].index_l[siguiente] = l;
    Bucket[dist].num_elem++;
    numelembucket[dist]++;
    if (maps[2][new_mapindex] > dreal) {
      maps[2][new_mapindex] = dreal;
    }
  }

  return 0;
}

int DToptimo3d(char *prototypes, int height, int width, int depth, int K, float **maps, int object_label, int domain_label, float hx, float hy, float hz) {
  int i, j, x, y, z, l, r, count;
  int siguiente, indice_actual;
  int d, mapindex, buckets_empty;
  struct element *Element;
  struct element *Proto;
  struct bucket Bucket[NUM_BUCKETS];
  float distancia_from_l, distancia_ultima;

  /* propagate3d hardcodes maps[2] as the real-distance accumulator, and
     assign() writes the integer bucket distance into maps[K]; those only
     stay distinct slots for K==1. Every caller also currently allocates
     maps with exactly 3 rows (0..2). Reject K>=2 here instead of silently
     corrupting maps[2] until both this function and its callers are
     generalized for K>1. */
  if (K >= 2) {
    printf("Error: DToptimo3d only supports K=1 (maps[2] collides with maps[K] for K>=2)\n");
    return 1;
  }

  /* inicializamos los mapas de etiquetas a -1 */
  for (i = 0; i < K; i++) {
    for (j = 0; j < height * width * depth; j++) {
      maps[i][j] = -1;
    }
  }
  for (j = 0; j < height * width * depth; j++) {
    maps[2][j] = 999999;
  }
  /* reservamos memoria para los prototipos, para los Elementos 
     y para el bucket inicial (bucket 0) */
  Proto = (struct element *)malloc(sizeof(struct element) * height * width * depth);
  Element = (struct element *)malloc(sizeof(struct element) * height * width * depth);
  Bucket[0].index_elem = (int *)malloc(sizeof(int) * MAX_ELEM_IN_BUCKET);
  Bucket[0].index_l = (int *)malloc(sizeof(int) * MAX_ELEM_IN_BUCKET);
  memset(Bucket[0].index_elem, 0, sizeof(int) * MAX_ELEM_IN_BUCKET);
  memset(Bucket[0].index_l, 0, sizeof(int) * MAX_ELEM_IN_BUCKET);
  Bucket[0].num_elem = 0;

  if (Element == (struct element *)NULL || Proto == (struct element *)NULL) {
    printf("Mucho me temo que te has quedado sin memoria \nCierra el Netscape, Kazaa, y todo lo que chupe memoria como una cochina,\ne intentalo otra vez\n");
    return 1;
  }

  /* inicializamos los elementos los prototipos, y el bucket inicial */
  i = 0;
  for (z = 0; z < depth; z++) {
    for (x = 0; x < height; x++) {
      for (y = 0; y < width; y++) {
        Element[i].icur = 0;
        Element[i].dcur = -1;
        Element[i].idcur = 0;
        Element[i].x = x;
        Element[i].y = y;
        Element[i].z = z;
        i++;
      }
    }
  }
  l = 0;
  count = 0;
  for (z = 0; z < depth; z++) {
    for (x = 0; x < height; x++) {
      for (y = 0; y < width; y++) {
        if (prototypes[count] == object_label) {
          Proto[l].x = x;
          Proto[l].y = y;
          Proto[l].z = z;
          Proto[l].icur = prototypes[count];
          Bucket[0].index_elem[l] = count;
          Bucket[0].index_l[l] = l;
          Bucket[0].num_elem++;
          l++;
        }
        count++;
      }
    }
  }

  /* printf("hay %d prototipos \n",l);*/
  for (i = 1; i < NUM_BUCKETS; i++) {
    Bucket[i].num_elem = 0;
    numelembucket[i] = 0;
  }
  numelembucket[0] = Bucket[0].num_elem;
  /* Fin de inicializacion */

  d = 0;
  while (1) {
    /*printf("Distancia actual: %d\n",d);*/
    if (d + 1 >= NUM_BUCKETS) {
      printf("excedido el maximo numero de buckets: %d\n", NUM_BUCKETS);
      return 1;
    }
    Bucket[d + 1].index_elem = (int *)malloc(sizeof(int) * MAX_ELEM_IN_BUCKET);
    Bucket[d + 1].index_l = (int *)malloc(sizeof(int) * MAX_ELEM_IN_BUCKET);
    /* printf("reservada memoria para bucket d+1\n" ); */
    while (Bucket[d].num_elem != 0) {
      /* Get (p,l) from Bucket d */
      /* printf("Obteniendo de bucket: %d\n",d); */
      siguiente = Bucket[d].num_elem - 1;
      /* printf("num_elem = %d \n",Bucket[d].num_elem); */
      indice_actual = Bucket[d].index_elem[siguiente];
      l = Bucket[d].index_l[siguiente];
      Bucket[d].index_elem[siguiente] = -1;
      Bucket[d].num_elem--;
      mapindex = mapIndex3D(Element[indice_actual].x, Element[indice_actual].y, Element[indice_actual].z, height, width, depth);
      if (mapindex != indice_actual) {
        printf("Error catastrofico: mapindex %d != indice_actual %d \n", mapindex, indice_actual);
        return 1;
      }
      /* printf("mapindex %d\n",mapindex);  */
      /*         *****         */

      if (Element[indice_actual].icur < K) {
        if (Element[indice_actual].dcur < d) {
          assign(&Element[indice_actual], l, mapindex, d, maps, K);
          /* printf("asignado pixel\n");*/
          if (propagate3d(Element[indice_actual], l, height, width, depth, Proto, Bucket, maps, d, Element, K, domain_label, prototypes, hx, hy, hz) != 0) return 1;
        } else {
          for (j = Element[indice_actual].idcur; j <= Element[indice_actual].icur; j++) {
            if (maps[j - 1][mapindex] == l) { break; }
            if (j == Element[indice_actual].icur) {
              distancia_from_l = distance3d_voxelsize(Element[indice_actual].x, Element[indice_actual].y, Element[indice_actual].z, Proto[l].x, Proto[l].y, Proto[l].z, hx, hy, hz);
              r = maps[j - 1][mapindex];
              distancia_ultima = distance3d_voxelsize(Element[indice_actual].x, Element[indice_actual].y, Element[indice_actual].z, Proto[r].x, Proto[r].y, Proto[r].z, hx, hy, hz);
              if (distancia_from_l < distancia_ultima) {
                maps[j - 1][mapindex] = l;
                assign(&Element[indice_actual], r, mapindex, d, maps, K);
              } else {
                assign(&Element[indice_actual], l, mapindex, d, maps, K);
              }
              if (propagate3d(Element[indice_actual], l, height, width, depth, Proto, Bucket, maps, d, Element, K, domain_label, prototypes, hx, hy, hz) != 0) return 1;
              break;
            }
          } /* end for */
        } /*end else*/
      } /* end if (Element_p.icur < K) */
    } /*end while */

    if (Bucket[d].num_elem == 0) {
      /* printf("num_elem en bucket %d = %d\n",d,numelembucket[d]);*/
      free(Bucket[d].index_elem);
      free(Bucket[d].index_l);
    }
    d++;
    if (d >= NUM_BUCKETS) {
      printf("excedido el maximo numero de buckets: %d\n", NUM_BUCKETS);
      return 1;
    }
    buckets_empty = 1;
    for (i = 0; i <= d; i++) {
      if (Bucket[d].num_elem != 0) {
        buckets_empty = 0;
        break;
      }
    }
    /* if all buckets are empty break the main bucle and finalize */
    if (buckets_empty == 1) break;
  }

  printf("Distancia maxima alcanzada: %d\n", d);

  /* liberamos memoria */
  free(Proto);
  free(Bucket[d].index_elem);
  free(Bucket[d].index_l);
  free(Element);
  printf("DT optimo3d OK\n");
  return 0;
}
