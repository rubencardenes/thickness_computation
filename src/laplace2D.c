/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "laplace2D.h"

/* Writes the solver's starting field from the label image: the band (label 2)
   starts at 0, the region outside the domain (label 255) at 255 -- or -1 when
   `reverse` flips which boundary the field runs towards -- and every other
   label is held at its own value as the Dirichlet condition the relaxation
   reads. Shared by laplace2D and poisson2D, which differ only in their update
   step. */
void init_laplace_field2D(const unsigned char *input, int height, int width, float **output, int reverse) {
  int i, j;
  int sum = 0;

  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      if (input[sum] == 2) {
        output[i][j] = 0;
      } else {
        if (input[sum] == 255) {
          output[i][j] = reverse ? -1 : 255;
        } else {
          output[i][j] = input[sum];
        }
      }
      sum++;
    }
  }
}

/* Input values
255 Outside domain
1 Exterior boundary 
0 Interior boundary
2 Inside domain 
*/
int laplace2D(unsigned char *input, int height, int width, float **output, int iterations, float lambda, int reverse) {
  int i, j, l;
  int sum = 0;
  init_laplace_field2D(input, height, width, output, reverse);

  /* Solve Laplacian */
  for (l = 0; l < iterations; l++) {
    sum = 0;
    for (i = 0; i < height; i++) {
      for (j = 0; j < width; j++) {
        if (input[sum] == 2 && i != 0 && i != height - 1 && j != 0 && j != width - 1) {
          output[i][j] = output[i][j] + (lambda + 1) * (0.25 * (output[i - 1][j] + output[i + 1][j] + output[i][j - 1] + output[i][j + 1]) - output[i][j]);
        }
        sum++;
      }
    }
  }

  return 0;
}

int EdgeDetect(unsigned char *domain, int height, int width) {
  int x, y, i;
  i = 0;

  for (x = 0; x < height; x++) {
    for (y = 0; y < width; y++) {
      if ((x == 0) || (y == 0) || (x == height - 1) || (y == width - 1)) {
        /* nothing to do */
      } else if ((domain[i] != 0) &&
                 ((domain[i + 1] == 0) ||
                  (domain[i - 1] == 0) ||

                  (domain[i + width] == 0) ||
                  (domain[i - width] == 0))) {

        domain[i] = 1;
      }
      i++;
    }
  }
  return 0;
}


/* Collects the whole 8-connected component of `label` that contains `seed`
   into `component` (cleared first), marking each collected pixel in `visited`.
   A pixel is marked when it is pushed, so it is collected exactly once and
   component->num_elem ends up as the component's size in pixels. The list is
   both the worklist and the record: list_take walks it without consuming it.

   This is the one definition of "component" in this file; maxcomponent2D and
   sizefilter2D must agree on it, since compute_boundary2D pairs them. */
static void collect_component2D(const unsigned short *data, int height, int width, int label,
                                int seed, unsigned char *visited, struct index_list *component) {
  int neighbors[MAX_NEIGHBORS_2D];
  int mapindex, k, count;

  list_clear(component);
  visited[seed] = 1;
  list_push(component, seed);
  while ((mapindex = list_take(component)) >= 0) {
    count = neighbors2D(mapindex, height, width, neighbors);
    for (k = 0; k < count; k++) {
      if (data[neighbors[k]] == label && !visited[neighbors[k]]) {
        visited[neighbors[k]] = 1;
        list_push(component, neighbors[k]);
      }
    }
  }
}

/* Returns the size in pixels of the largest 8-connected component of value
   `label` in `data`, 0 if no pixel carries `label`, or -1 if it could not
   allocate. */
int maxcomponent2D(unsigned short *data, int height, int width, int label) {
  struct index_list component;
  int npixels = height * width;
  int seed, max_size;
  unsigned char *visited;

  visited = (unsigned char *)calloc(npixels, sizeof(unsigned char));
  if (visited == NULL) {
    fprintf(stderr, "maxcomponent2D: out of memory\n");
    return -1;
  }
  if (list_init(&component, npixels) != 0) {
    fprintf(stderr, "maxcomponent2D: out of memory\n");
    free(visited);
    return -1;
  }

  max_size = 0;
  for (seed = 0; seed < npixels; seed++) {
    if (data[seed] != label || visited[seed]) continue;
    collect_component2D(data, height, width, label, seed, visited, &component);
    if (component.num_elem > max_size) {
      max_size = component.num_elem;
    }
  }

  list_free(&component);
  free(visited);

  return max_size;
}

/* Relabels every 8-connected component of `oldlabel` whose size is `max_size`
   pixels or fewer to `newlabel`; larger components keep `oldlabel`. Returns 0,
   or 1 if it could not allocate. */
int sizefilter2D(unsigned short *data, int height, int width, int max_size, int oldlabel, int newlabel) {
  struct index_list component;
  int npixels = height * width;
  int seed, k;
  unsigned char *visited;

  visited = (unsigned char *)calloc(npixels, sizeof(unsigned char));
  if (visited == NULL) {
    fprintf(stderr, "sizefilter2D: out of memory\n");
    return 1;
  }
  if (list_init(&component, npixels) != 0) {
    fprintf(stderr, "sizefilter2D: out of memory\n");
    free(visited);
    return 1;
  }

  for (seed = 0; seed < npixels; seed++) {
    if (data[seed] != oldlabel || visited[seed]) continue;
    collect_component2D(data, height, width, oldlabel, seed, visited, &component);
    if (component.num_elem <= max_size) {
      for (k = 0; k < component.num_elem; k++) {
        data[component.elem[k]] = newlabel;
      }
    }
  }

  list_free(&component);
  free(visited);

  return 0;
}

/* Relabels the 8-connected region of `oldlabel` reachable from `startindex` to
   `newlabel`. Returns 0, or 1 if it could not allocate.

   Iterative rather than recursive: the recursive version descended once per
   pixel, so a large region could exhaust the call stack. A pixel is relabeled
   as it is pushed, which is also what stops it from being pushed twice. */
int floodfill(unsigned char *domain, int startindex, unsigned short oldlabel, unsigned char newlabel, int height, int width) {
  int neighbors[MAX_NEIGHBORS_2D];
  struct index_list region;
  int mapindex, k, count;

  if (domain[startindex] != oldlabel || oldlabel == newlabel) {
    return 0;
  }
  if (list_init(&region, height * width) != 0) {
    fprintf(stderr, "floodfill: out of memory\n");
    return 1;
  }

  domain[startindex] = newlabel;
  list_push(&region, startindex);
  while ((mapindex = list_take(&region)) >= 0) {
    count = neighbors2D(mapindex, height, width, neighbors);
    for (k = 0; k < count; k++) {
      if (domain[neighbors[k]] == oldlabel) {
        domain[neighbors[k]] = newlabel;
        list_push(&region, neighbors[k]);
      }
    }
  }

  list_free(&region);
  return 0;
}

/* Splits the two boundaries of the domain into the Dirichlet values the Laplace
   solver reads: the first 8-connected run of label-1 pixels becomes 0, and the
   region reachable from pixel 0 (the background outside the domain) becomes 3.
   Returns 0, or 1 if no boundary pixel was found or it could not allocate. */
int RelabelBoundary(unsigned char *domain, int height, int width) {
  int neighbors[MAX_NEIGHBORS_2D];
  struct index_list boundary;
  int npixels = height * width;
  int mapindex, k, count, start;

  start = -1;
  for (mapindex = 0; mapindex < npixels; mapindex++) {
    if (domain[mapindex] == 1) {
      start = mapindex;
      break;
    }
  }

  /* No boundary pixel (value 1) was found: the input domain does not have the
     expected form (a non-zero region bordering a zero background), so there is
     nothing to relabel. Bail out instead of writing to an invalid index. */
  if (start < 0) {
    fprintf(stderr, "RelabelBoundary: no boundary pixel found; check input labels\n");
    return 1;
  }

  if (list_init(&boundary, npixels) != 0) {
    fprintf(stderr, "RelabelBoundary: out of memory\n");
    return 1;
  }

  domain[start] = 0;
  list_push(&boundary, start);
  while ((mapindex = list_take(&boundary)) >= 0) {
    count = neighbors2D(mapindex, height, width, neighbors);
    for (k = 0; k < count; k++) {
      if (domain[neighbors[k]] == 1) {
        domain[neighbors[k]] = 0;
        list_push(&boundary, neighbors[k]);
      }
    }
  }
  list_free(&boundary);

  /* floodfill the exterior region, relabeling to 3 */
  return floodfill(domain, 0, domain[0], 3, height, width);
}


int iGradY(float **ppfData, float **ppfGradient, int numRowX, int numColY)
/* 
 * PURPOSE : Calculate the partial derivative in X direction.
 * 			X - increases across the rows.
 * 			Y - increases across the columns.
 * INPUTS :
 * ppfData     : The image data.
 * ppfGradient : The gradient will be stored in this array
 * numRowX	   : The number of points in the X direction.
 * numColY	   : The number of points in the Y direction.
 *
 * OUTPUTS :
 * 		Ix is approximated.
 *
 * METHOD : At each point in the image a central difference approx.
 * 	    to the first derivative in the x direction is made.
 * 	    At the boundaries a forward difference or a backward difference
 * 	    approximation is made.
 */
{
  register int i;
  register int j;
  register int iMaxX;
  register int iMaxY;

  iMaxX = numRowX - 1;
  iMaxY = numColY - 1;
  assert(iMaxX >= 0 && iMaxY >= 0);

  if (iMaxX == 0) { /* Ix == 0 everywhere for one row data */
    for (j = 0; j <= iMaxY; j++) {
      ppfGradient[0][j] = 0;
    }
  } else { /* iMaxY can still be zero - but this is OK */
    // Along edges of image where i == 0 and i == iMaxX (top and bottom row)
    //						       FDX			BDX
    // Since a central difference Ix approximation is not possible.
    // This reaches all four corners.
    for (j = 0; j <= iMaxY; j++) {
      ppfGradient[0][j] = FDX(ppfData, 0, j);
      ppfGradient[iMaxX][j] = BDX(ppfData, iMaxX, j);
    }

    // On all internal pixels we can calculate a central difference approx.
    for (i = 1; i < iMaxX; i++) {
      for (j = 1; j < iMaxY; j++) {
        ppfGradient[i][j] = CDX(ppfData, i, j);
      }
    }
    // Along edges of image where j == 0 and j == iMaxY (left, right col).
    // we can calculate CDX, but not in the corners.
    for (i = 1; i < iMaxX; i++) {
      ppfGradient[i][0] = CDX(ppfData, i, 0);
      ppfGradient[i][iMaxY] = CDX(ppfData, i, iMaxY);
    }
  }
  return 0; // Success
}


int iGradX(float **ppfData, float **ppfGradient, int numRowX, int numColY)
/* 
 * PURPOSE : Calculate the partial derivative in Y direction.
 * 			X - increases across the rows.
 * 			Y - increases across the columns.
 * INPUTS :
 * ppfData     : The image data.
 * ppfGradient : The gradient will be stored in this array
 * numRowX	   : The number of points in the X direction.
 * numColY	   : The number of points in the Y direction.
 *
 * OUTPUTS :
 * 		Iy is approximated.
 *
 * METHOD : At each point in the image a central difference approx.
 * 	    to the first derivative in the y direction is made.
 * 	    At the boundaries a forward difference or a backward difference
 * 	    approximation is made.
 */
{
  register int i;
  register int j;
  register int iMaxX;
  register int iMaxY;

  iMaxX = numRowX - 1;
  iMaxY = numColY - 1;
  assert(iMaxX >= 0 && iMaxY >= 0);

  if (iMaxY == 0) { /* Iy == 0 everywhere for one column data */
    for (i = 0; i <= iMaxX; i++) {
      ppfGradient[i][0] = 0;
    }
  } else { /* iMaxX can still be zero - but this is OK */
    // Along edges of image where j == 0 and j == iMaxY (left and right col)
    //						       FDY			BDY
    // Since a central difference Iy approximation is not possible.
    // This reaches all four corners.
    for (i = 0; i <= iMaxX; i++) {
      ppfGradient[i][0] = FDY(ppfData, i, 0);
      ppfGradient[i][iMaxY] = BDY(ppfData, i, iMaxY);
    }

    // On all internal pixels we can calculate a central difference approx.
    for (i = 1; i < iMaxX; i++) {
      for (j = 1; j < iMaxY; j++) {
        ppfGradient[i][j] = CDY(ppfData, i, j);
      }
    }
    // Along edges of image where i == 0 and i == iMaxX (top, bottom rows).
    // we can calculate CDY, but not in the corners.
    for (j = 1; j < iMaxY; j++) {
      ppfGradient[0][j] = CDY(ppfData, 0, j);
      ppfGradient[iMaxX][j] = CDY(ppfData, iMaxX, j);
    }
  }
  return 0; // Success
}

int normalize(float **gradientx, float **gradienty, int height, int width) {
  int i, j;
  float norma;

  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      norma = sqrt(gradientx[i][j] * gradientx[i][j] + gradienty[i][j] * gradienty[i][j]);
      if (norma != 0) {
        gradientx[i][j] = gradientx[i][j] / norma;
        gradienty[i][j] = gradienty[i][j] / norma;
      }
    }
  }
  return 0;
}

/* Upper bound on the corners this detects. `tramo` counts the arcs the corners
   cut the curve into and is used to index count0/count1, so those need one more
   slot than there are corners. */
#define MAX_CORNERS 10

int new_compute_corners(unsigned short *input, int height, int width) {
  int i, k, start, end, mapindex, newmapindex, flag, num_corners = 0;
  int in, ip, tramo, corner[MAX_CORNERS];
  int count, count0[MAX_CORNERS + 1], count1[MAX_CORNERS + 1];
  float ax, ay, bx, by;
  unsigned char *aux;
  int neighbors[MAX_NEIGHBORS_2D];
  struct index_list curve;    /* the boundary walk, in traversal order */
  struct index_list frontier; /* pixels whose neighbours are still to be visited */
  int max_num_list = 2000, num_vec = 10;
  int *ideal;
  float *coseno_ideal;
  float coseno[num_vec + 1];

  aux = (unsigned char *)calloc(height * width, sizeof(unsigned char));
  if (aux == NULL || list_init(&curve, max_num_list) != 0 || list_init(&frontier, max_num_list) != 0) {
    fprintf(stderr, "new_compute_corners: out of memory\n");
    return 1;
  }

  i = 0;
  while (input[i] != 1 && i < height * width) {
    i++;
  }
  start = i;
  if (start == height * width) {
    /* printf("desapareciendo del mapa \n");*/
    list_free(&curve);
    list_free(&frontier);
    free(aux);
    return 0;
  }
  aux[start] = 1;
  list_push(&curve, start);
  list_push(&frontier, start);

  /* walk the label-1 side of the boundary, in FIFO order so that `curve` comes
     out ordered along the curve -- the curvature pass below reads it that way */
  while ((mapindex = list_take(&frontier)) >= 0) {
    count = neighbors2D(mapindex, height, width, neighbors);
    for (k = 0; k < count; k++) {
      newmapindex = neighbors[k];
      if (input[newmapindex] == 1 && aux[newmapindex] == 0) {
        aux[newmapindex] = 1;
        if (list_push(&curve, newmapindex) != 0 || list_push(&frontier, newmapindex) != 0) {
          printf("Error curve.num_elem > %d \n", max_num_list);
          list_free(&curve);
          list_free(&frontier);
          free(aux);
          return 1;
        }
      }
    }
  }
  end = curve.elem[curve.num_elem - 1];

  /* printf("num_elem en la curva %d\n",curve.num_elem);*/
  /* printf("start x %d y %d\n",maptox(start,width),maptoy(start,width));
  printf("end x %d y %d\n",maptox(end,width),maptoy(end,width));
  printf("input[end] %d input[start] %d\n",input[end],input[start]);*/

  /* tengo que crear una lista de puntos pertenecientes a la curva
     correspondientes a la otra parte*/
  list_clear(&frontier);
  list_push(&frontier, end);
  while ((mapindex = list_take(&frontier)) >= 0) {
    count = neighbors2D(mapindex, height, width, neighbors);
    for (k = 0; k < count; k++) {
      newmapindex = neighbors[k];
      if (input[newmapindex] == 0 && aux[newmapindex] == 0) {
        aux[newmapindex] = 1;
        if (list_push(&curve, newmapindex) != 0 || list_push(&frontier, newmapindex) != 0) {
          printf("Error curve.num_elem > %d \n", max_num_list);
          list_free(&curve);
          list_free(&frontier);
          free(aux);
          return 1;
        }
      }
    }
  }

  printf("num_elem en la curva %d\n", curve.num_elem);
  coseno_ideal = (float *)malloc(sizeof(float) * curve.num_elem);
  ideal = (int *)malloc(sizeof(int) * curve.num_elem);
  /* tengo que calcular (ax,ay), (bx,by) en una vecindad de num_vec puntos*/
  for (i = 0; i < curve.num_elem; i++) {
    for (k = num_vec; k > 0; k--) {
      ip = (i + k) % curve.num_elem;
      in = (i - k + curve.num_elem) % curve.num_elem;
      ax = maptox(curve.elem[i], width) - maptox(curve.elem[ip], width);
      ay = maptoy(curve.elem[i], width) - maptoy(curve.elem[ip], width);
      bx = maptox(curve.elem[i], width) - maptox(curve.elem[in], width);
      by = maptoy(curve.elem[i], width) - maptoy(curve.elem[in], width);
      coseno[k] = (float)(ax * bx + ay * by) / (sqrt(ax * ax + ay * ay) * sqrt(bx * bx + by * by));
      if (k < num_vec) {
        if (coseno[k] >= coseno[k + 1]) {
          ideal[i] = k + 1;
          coseno_ideal[i] = (float)coseno[k + 1];
          k = 0;
        }
      } else {
        ideal[i] = k;
        coseno_ideal[i] = coseno[k];
      }
    }
  }
  /* Necesito otra pasada para encontrar los puntos de max curvatura*/
  for (i = 0; i < curve.num_elem; i++) {
    flag = 1;
    for (k = -num_vec / 2; k <= num_vec / 2; k++) {
      if (k == 0) continue;
      in = (i + k + curve.num_elem) % curve.num_elem;
      if (coseno_ideal[in] >= coseno_ideal[i]) {
        flag = 0;
      }
    }
    if (flag == 1 && coseno_ideal[i] > 0) {
      printf("found corner at x %d y %d coseno_ideal[i] %f\n", maptox(curve.elem[i], width), maptoy(curve.elem[i], width), coseno_ideal[i]);
      /* Check before storing: the guard used to run after the write, so the
         (MAX_CORNERS+1)-th corner wrote one past the end of corner[]. */
      if (num_corners >= MAX_CORNERS) {
        printf("Error, more than %d corners found\n", MAX_CORNERS);
        free(coseno_ideal);
        free(ideal);
        free(aux);
        list_free(&curve);
        list_free(&frontier);
        return 1;
      }
      aux[curve.elem[i]] = 255;
      corner[num_corners] = i;
      num_corners++;
      input[curve.elem[i]] = 7;
    }
    /* input[curve.elem[i]] = i;*/
  }

  if (num_corners < 2) {
    free(coseno_ideal);
    free(ideal);
    free(aux);
    list_free(&curve);
    list_free(&frontier);
    return 0;
  }

  for (i = 0; i <= MAX_CORNERS; i++) {
    count0[i] = 0;
    count1[i] = 0;
  }
  k = (corner[0] + 1) % curve.num_elem;
  tramo = 0;
  for (i = 0; i < curve.num_elem; i++) {
    mapindex = curve.elem[k];
    if (aux[mapindex] == 255) {
      tramo++;
    }
    if (input[mapindex] == 0) {
      count0[tramo]++;
    }
    if (input[mapindex] == 1) {
      count1[tramo]++;
    }
    k = (k + 1) % curve.num_elem;
  }
  /* printf("count0[0] %d count1[0] %d\n",count0[0],count1[0]);
     printf("count0[1] %d count1[1] %d\n",count0[1],count1[1]);*/

  if (num_corners == 2) {
    if (count0[0] > count0[1]) {
      count0[0] = 0;
      count0[1] = 1;
    } else {
      count0[0] = 1;
      count0[1] = 0;
    }
    k = (corner[0] + 1) % curve.num_elem;
    tramo = 0;
    for (i = 0; i < curve.num_elem; i++) {
      mapindex = curve.elem[k];
      if (aux[mapindex] == 255) {
        tramo++;
      }
      input[mapindex] = count0[tramo];
      k = (k + 1) % curve.num_elem;
    }
  } else {
    k = (corner[0] + 1) % curve.num_elem;
    tramo = 0;
    for (i = 0; i < curve.num_elem; i++) {
      mapindex = curve.elem[k];
      if (aux[mapindex] == 255) {
        tramo++;
      }
      if (count0[tramo] > count1[tramo]) {
        input[mapindex] = 0;
      }
      if (count0[tramo] < count1[tramo]) {
        input[mapindex] = 1;
      }
      k = (k + 1) % curve.num_elem;
    }
  }

  free(coseno_ideal);
  free(ideal);
  free(aux);
  list_free(&curve);
  list_free(&frontier);

  return 0;
}
