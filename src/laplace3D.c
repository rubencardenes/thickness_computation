/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "laplace3D.h"

/* Input values
255 Outside domain
1 Exterior boundary 
0 Interior boundary
2 Inside domain 
*/
int laplace3D(unsigned char *input, int height, int width, int depth, float ***output, int iterations, float lambda) {
  int i, j, k, l;
  int sum = 0;
  /* Initialize domain, inside=0, and boundaries values.
     Iteration order (k,j,i with i innermost) must match the i-fastest/
     j-next(height stride)/k-slowest(height*width stride) convention that
     compute_boundary_cortex3D/EdgeDetect3D use to populate `input`. */
  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        if (input[sum] == 2) {
          output[k][i][j] = 0;
        } else if (input[sum] == 3) {
          output[k][i][j] = -1;
        } else {
          output[k][i][j] = input[sum];
        }
        sum++;
      }
    }
  }

  /* Solve Laplacian */
  for (l = 0; l < iterations; l++) {
    sum = 0;
    for (k = 0; k < depth; k++) {
      for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++) {
          if (input[sum] == 2 && i != 0 && i != height - 1 && j != 0 && j != width - 1 && k != 0 && k != depth - 1) {
            /*output[k][i][j] = (output[k-1][i][j] + output[k+1][i][j] + output[k][i-1][j] + output[k][i+1][j] + output[k][i][j-1] + output[k][i][j+1])/6;*/
            output[k][i][j] = output[k][i][j] + (lambda + 1) * (0.16666667 * (output[k - 1][i][j] + output[k + 1][i][j] + output[k][i - 1][j] + output[k][i + 1][j] + output[k][i][j - 1] + output[k][i][j + 1]) - output[k][i][j]);
          }
          sum++;
        }
      }
    }
  }

  return 0;
}

int laplace3D_voxelsize(unsigned char *input, int height, int width, int depth, float ***output, int iterations, float hx, float hy, float hz, float lambda) {
  int i, j, k, l;
  int sum = 0;
  /* Initialize domain, inside=0, and boundaries values.
     Iteration order (k,j,i with i innermost) must match the i-fastest/
     j-next(height stride)/k-slowest(height*width stride) convention that
     compute_boundary_cortex3D/EdgeDetect3D use to populate `input`. */
  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        if (input[sum] == 2) {
          output[k][i][j] = 0;
        } else if (input[sum] == 3) {
          output[k][i][j] = -1;
        } else {
          output[k][i][j] = input[sum];
        }
        sum++;
      }
    }
  }

  /* Solve Laplacian */
  for (l = 0; l < iterations; l++) {
    sum = 0;
    for (k = 0; k < depth; k++) {
      for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++) {
          if (input[sum] == 2 && i != 0 && i != height - 1 && j != 0 && j != width - 1 && k != 0 && k != depth - 1) {
            output[k][i][j] = 0.5 * ((output[k - 1][i][j] + output[k + 1][i][j]) / (hz * hz) + (output[k][i - 1][j] + output[k][i + 1][j]) / (hy * hy) + (output[k][i][j - 1] + output[k][i][j + 1]) / (hx * hx)) * (hx * hx * hy * hy * hz * hz) / (hx * hx * hy * hy + hy * hy * hz * hz + hx * hx * hz * hz);
            /* output[k][i][j] = output[k][i][j]+(lambda+1)*(0.5 *( (output[k-1][i][j] + output[k+1][i][j])/(hz*hz) + (output[k][i-1][j] + output[k][i+1][j])/(hy*hy) + (output[k][i][j-1] + output[k][i][j+1])/(hx*hx))*(hx*hx*hy*hy*hz*hz)/(hx*hx*hy*hy + hy*hy*hz*hz + hx*hx*hz*hz) - output[k][i][j]);*/
          }
          sum++;
        }
      }
    }
  }

  return 0;
}

int EdgeDetect3D_knee(unsigned char *domain, int height, int width, int depth) {
  int x, y, z, i;
  i = 0;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < width; y++) {
      for (x = 0; x < height; x++) {
        if ((x == 0) || (y == 0) || (z == 0) || (x == height - 1) || (y == width - 1) || (z == depth - 1)) {
          /* domain[i]=255;*/
        } else if ((domain[i] == 2) &&
                   ((domain[i + 1] == 255) ||
                    (domain[i - 1] == 255) ||

                    (domain[i + height] == 255) ||
                    (domain[i - height] == 255) ||

                    (domain[i + height * width] == 255) ||
                    (domain[i - height * width] == 255))) {

          domain[i] = 1;
        }
        /*else {	   	   
    domain[i]=255;
    }*/
        i++;
      }
    }
  }
  return 0;
}

int EdgeDetect3D(unsigned char *domain, int height, int width, int depth) {
  int x, y, z, i;
  i = 0;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < width; y++) {
      for (x = 0; x < height; x++) {
        if ((x == 0) || (y == 0) || (z == 0) || (x == height - 1) || (y == width - 1) || (z == depth - 1)) {
          /* domain[i]=255;*/
        } else if ((domain[i] != 0) &&
                   ((domain[i + 1] == 0) ||
                    (domain[i - 1] == 0) ||

                    (domain[i + height] == 0) ||
                    (domain[i - height] == 0) ||

                    (domain[i + height * width] == 0) ||
                    (domain[i - height * width] == 0))) {

          domain[i] = 1;
        }
        /*else {	   	   
    domain[i]=255;
    }*/
        i++;
      }
    }
  }
  return 0;
}

/* Relabels the 6-connected region of `oldlabel` reachable from `startindex` to
   `newlabel` (face neighbours only, unlike the 26-connected RelabelBoundary3D).
   Returns 0, or 1 if it could not allocate. A voxel is relabeled as it is
   pushed, which is what stops it from being pushed twice. */
int floodfill3D(unsigned char *domain, int startindex, unsigned short oldlabel, unsigned char newlabel, int height, int width, int depth) {
  int neighbors[MAX_NEIGHBORS_3D_FACES];
  struct index_list region;
  int mapindex, k, count;

  if (domain[startindex] != oldlabel || oldlabel == newlabel) {
    return 0;
  }
  if (list_init(&region, height * width * depth) != 0) {
    fprintf(stderr, "floodfill3D: out of memory\n");
    return 1;
  }

  domain[startindex] = newlabel;
  list_push(&region, startindex);
  while ((mapindex = list_take(&region)) >= 0) {
    count = neighbors3D_faces(mapindex, height, width, depth, neighbors);
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

/* 3D counterpart of RelabelBoundary: the first 26-connected run of label-1
   voxels becomes 0, then the region reachable from voxel 0 (the background
   outside the domain) becomes 3. Returns 0, or 1 if no boundary voxel was
   found or it could not allocate. */
int RelabelBoundary3D(unsigned char *domain, int height, int width, int depth) {
  int neighbors[MAX_NEIGHBORS_3D];
  struct index_list boundary;
  int nvoxels = height * width * depth;
  int mapindex, k, count, start;

  start = -1;
  for (mapindex = 0; mapindex < nvoxels; mapindex++) {
    if (domain[mapindex] == 1) {
      start = mapindex;
      break;
    }
  }

  if (start < 0) {
    fprintf(stderr, "RelabelBoundary3D: no boundary voxel found; check input labels\n");
    return 1;
  }

  if (list_init(&boundary, nvoxels) != 0) {
    fprintf(stderr, "RelabelBoundary3D: out of memory\n");
    return 1;
  }

  domain[start] = 0;
  list_push(&boundary, start);
  while ((mapindex = list_take(&boundary)) >= 0) {
    count = neighbors3D(mapindex, height, width, depth, neighbors);
    for (k = 0; k < count; k++) {
      if (domain[neighbors[k]] == 1) {
        domain[neighbors[k]] = 0;
        list_push(&boundary, neighbors[k]);
      }
    }
  }
  list_free(&boundary);

  printf("Doing floodfill3D, domain[0] %d \n", domain[0]);
  /* floodfill the exterior region, relabeling to 3 */
  return floodfill3D(domain, 0, domain[0], 3, height, width, depth);
}


int iGradX3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY, int numSlice, float hx)
/* 
 * PURPOSE : Calculate the partial derivative in X direction.
 * 			X - increases across the rows.
 * 			Y - increases across the columns.
 * 			Z - increases across the slices.
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
  register int k;
  register int iMaxX;
  register int iMaxY;
  register int iMaxZ;

  iMaxX = numRowX - 1;
  iMaxY = numColY - 1;
  iMaxZ = numSlice - 1;
  assert(iMaxX >= 0 && iMaxY >= 0 && iMaxZ >= 0);

  if (iMaxX == 0) { /* Ix == 0 everywhere for one row data */
    for (k = 0; k <= iMaxZ; k++) {
      for (i = 0; i <= iMaxY; i++) {
        ppfGradient[k][i][0] = 0;
      }
    }
  } else { /* iMaxY can still be zero - but this is OK */
           // Along edges of image where i == 0 and i == iMaxX (top and bottom row)
           //						       FDX			BDX
           // Since a central difference Ix approximation is not possible.
           // This reaches all four corners.
    for (k = 0; k <= iMaxZ; k++) {
      for (i = 0; i <= iMaxY; i++) {
        ppfGradient[k][i][0] = FDX3(ppfData, k, i, 0, hx);
        ppfGradient[k][i][iMaxX] = BDX3(ppfData, k, i, iMaxX, hx);
      }
    }

    // On all internal pixels we can calculate a central difference approx.
    for (k = 0; k <= iMaxZ; k++) {
      for (i = 0; i <= iMaxY; i++) {
        for (j = 1; j < iMaxX; j++) {
          ppfGradient[k][i][j] = CDX3(ppfData, k, i, j, hx);
        }
      }
    }
  }
  return 0; // Success
}


int iGradY3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY, int numSlice, float hy)
/* 
 * PURPOSE : Calculate the partial derivative in Y direction.
 * 			X - increases across the columns.
 * 			Y - increases across the rows. 
 * 			Z - increases across the slices.
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
  register int k;
  register int iMaxX;
  register int iMaxY;
  register int iMaxZ;

  iMaxX = numRowX - 1;
  iMaxY = numColY - 1;
  iMaxZ = numSlice - 1;
  assert(iMaxX >= 0 && iMaxY >= 0 && iMaxZ >= 0);

  if (iMaxY == 0) { /* Iy == 0 everywhere for one column data */
    for (k = 0; k <= iMaxZ; k++) {
      for (j = 0; j <= iMaxX; j++) {
        ppfGradient[k][0][j] = 0;
      }
    }
  } else { /* iMaxX can still be zero - but this is OK */
           // Along edges of image where j == 0 and j == iMaxY (left and right col)
           //						       FDY			BDY
           // Since a central difference Iy approximation is not possible.
           // This reaches all four corners.
    for (k = 0; k <= iMaxZ; k++) {
      for (j = 0; j <= iMaxX; j++) {
        ppfGradient[k][0][j] = FDY3(ppfData, k, 0, j, hy);
        ppfGradient[k][iMaxY][j] = BDY3(ppfData, k, iMaxY, j, hy);
      }
    }

    // On all internal pixels we can calculate a central difference approx.
    for (k = 0; k <= iMaxZ; k++) {
      for (i = 1; i < iMaxY; i++) {
        for (j = 0; j <= iMaxX; j++) {
          ppfGradient[k][i][j] = CDY3(ppfData, k, i, j, hy);
        }
      }
    }
  }
  return 0; // Success
}

int iGradZ3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY, int numSlice, float hz)
/* 
 * PURPOSE : Calculate the partial derivative in Z direction.
 * 			X - increases across the columns.
 * 			Y - increases across the rows. 
 * 			Z - increases across the slices.
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
  register int k;
  register int iMaxX;
  register int iMaxY;
  register int iMaxZ;

  iMaxX = numRowX - 1;
  iMaxY = numColY - 1;
  iMaxZ = numSlice - 1;
  assert(iMaxX >= 0 && iMaxY >= 0 && iMaxZ >= 0);

  if (iMaxZ == 0) { /* Iy == 0 everywhere for one column data */
    for (i = 0; i <= iMaxY; i++) {
      for (j = 0; j <= iMaxX; j++) {
        ppfGradient[0][i][j] = 0;
      }
    }
  } else { /* iMaxX can still be zero - but this is OK */
           // Along edges of image where j == 0 and j == iMaxY (left and right col)
           //						       FDY			BDY
           // Since a central difference Iy approximation is not possible.
           // This reaches all four corners.
    for (i = 0; i <= iMaxY; i++) {
      for (j = 0; j <= iMaxX; j++) {
        ppfGradient[0][i][j] = FDZ3(ppfData, 0, i, j, hz);
        ppfGradient[iMaxZ][i][j] = BDZ3(ppfData, iMaxZ, i, j, hz);
      }
    }

    // On all internal pixels we can calculate a central difference approx.
    for (k = 1; k < iMaxZ; k++) {
      for (i = 0; i <= iMaxY; i++) {
        for (j = 0; j <= iMaxX; j++) {
          ppfGradient[k][i][j] = CDZ3(ppfData, k, i, j, hz);
        }
      }
    }
  }
  return 0; // Success
}

int normalize3D(float ***gradientx, float ***gradienty, float ***gradientz, int height, int width, int depth) {
  int i, j, k;
  float norma;

  for (k = 0; k < depth; k++) {
    for (i = 0; i < height; i++) {
      for (j = 0; j < width; j++) {
        norma = sqrt(gradientx[k][i][j] * gradientx[k][i][j] + gradienty[k][i][j] * gradienty[k][i][j] + gradientz[k][i][j] * gradientz[k][i][j]);
        if (norma != 0) {
          gradientx[k][i][j] = gradientx[k][i][j] / norma;
          gradienty[k][i][j] = gradienty[k][i][j] / norma;
          gradientz[k][i][j] = gradientz[k][i][j] / norma;
        }
      }
    }
  }
  return 0;
}
