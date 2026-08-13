/* Copyright (c) Ruben Cardenes Almeida 22/03/2002 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "thickness3D.h"
#include "laplace3D.h"

#define INF 9999999

/* 3D Yezzi PDE thickness propagation, reverse direction, "relaxed" variant:
   combines whichever of the up to 6 axis-neighbor values (x/y/z) are already
   propagated, with no per-axis tolerance rejection (unlike distanceYezzi_reverse3D).
   Used on the first pass, when few neighbors are populated yet. Returns -1 if
   none of the 6 neighbors have been reached.

   An axis whose gradient component is exactly 0 is skipped entirely: it means
   no transport along that axis, so its neighbours carry no information. It
   contributes fabs(0)*value == 0 to the numerator and 0 to the denominator, so
   letting it set `flag` would accept a voxel on the strength of an axis that
   said nothing -- and the value then computed comes from no upstream data at
   all. The forward variant used to do exactly that (a two-way `else` instead of
   `else if (g > 0)`); both now skip zero. */
float distanceYezzi_reverse3D_relax(float ***gradientx, float ***gradienty, float ***gradientz, int newmapindex, int x, int y, int z, float *maps, int height, int width, float r, float hx, float hy, float hz) {
  float distf;
  int flag = 0;

  distf = 1.0;
  if (gradientx[z][y][x] > 0) {
    if (maps[newmapindex + 1] > -1) {
      distf += fabs(gradientx[z][y][x]) * maps[newmapindex + 1] / hx;
      flag = 1;
    }
  } else if (gradientx[z][y][x] < 0) {
    if (maps[newmapindex - 1] > -1) {
      distf += fabs(gradientx[z][y][x]) * maps[newmapindex - 1] / hx;
      flag = 1;
    }
  }

  if (gradienty[z][y][x] > 0) {
    if (maps[newmapindex + width] > -1) {
      distf += fabs(gradienty[z][y][x]) * maps[newmapindex + width] / hy;
      flag = 1;
    }
  } else if (gradienty[z][y][x] < 0) {
    if (maps[newmapindex - width] > -1) {
      distf += fabs(gradienty[z][y][x]) * maps[newmapindex - width] / hy;
      flag = 1;
    }
  }

  if (gradientz[z][y][x] > 0) {
    if (maps[newmapindex + height * width] > -1) {
      distf += fabs(gradientz[z][y][x]) * maps[newmapindex + height * width] / hz;
      flag = 1;
    }
  } else if (gradientz[z][y][x] < 0) {
    if (maps[newmapindex - height * width] > -1) {
      distf += fabs(gradientz[z][y][x]) * maps[newmapindex - height * width] / hz;
      flag = 1;
    }
  }

  if (flag == 0) {
    distf = -1;
  } else {
    distf = distf / (fabs(gradientx[z][y][x]) / hx + fabs(gradienty[z][y][x]) / hy + fabs(gradientz[z][y][x]) / hz);
  }

  return distf;
}

/* Forward-direction counterpart of distanceYezzi_reverse3D_relax, including the
   same zero-gradient rule. */
float distanceYezzi3D_relax(float ***gradientx, float ***gradienty, float ***gradientz, int newmapindex, int x, int y, int z, float *maps, int height, int width, float r, float hx, float hy, float hz) {
  float distf;
  int flag = 0;

  distf = 1.0;
  if (gradientx[z][y][x] < 0) {
    if (maps[newmapindex + 1] > -1) {
      flag = 1;
      distf += fabs(gradientx[z][y][x]) * maps[newmapindex + 1] / hx;
    }
  } else if (gradientx[z][y][x] > 0) {
    if (maps[newmapindex - 1] > -1) {
      flag = 1;
      distf += fabs(gradientx[z][y][x]) * maps[newmapindex - 1] / hx;
    }
  }
  if (gradienty[z][y][x] < 0) {
    if (maps[newmapindex + width] > -1) {
      distf += fabs(gradienty[z][y][x]) * maps[newmapindex + width] / hy;
      flag = 1;
    }
  } else if (gradienty[z][y][x] > 0) {
    if (maps[newmapindex - width] > -1) {
      distf += fabs(gradienty[z][y][x]) * maps[newmapindex - width] / hy;
      flag = 1;
    }
  }

  if (gradientz[z][y][x] < 0) {
    if (maps[newmapindex + height * width] > -1) {
      distf += fabs(gradientz[z][y][x]) * maps[newmapindex + height * width] / hz;
      flag = 1;
    }
  } else if (gradientz[z][y][x] > 0) {
    if (maps[newmapindex - height * width] > -1) {
      distf += fabs(gradientz[z][y][x]) * maps[newmapindex - height * width] / hz;
      flag = 1;
    }
  }

  if (flag == 0) {
    distf = -1;
  } else {
    distf = distf / (fabs(gradientx[z][y][x]) / hx + fabs(gradienty[z][y][x]) / hy + fabs(gradientz[z][y][x]) / hz);
  }
  return distf;
}

/* Estimates the thickness value at (x,y,z) from the upstream x/y/z neighbors
   indicated by the Laplace gradient, rejecting the update (-1) if an upstream
   neighbor hasn't been reached yet and its gradient component exceeds tolerance
   r. Used by thickness3DYezzi_reverse once the front is past its first, relaxed
   pass. Unlike the 2D distanceYezzi_reverse, this one really does apply `r`. */
float distanceYezzi_reverse3D(float ***gradientx, float ***gradienty, float ***gradientz, int newmapindex, int x, int y, int z, float *maps, int height, int width, float r, float hx, float hy, float hz) {
  float distf;

  if (gradientx[z][y][x] > 0) {
    if (maps[newmapindex + 1] == -1 && fabs(gradientx[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex + 1] == -1) {
        distf = 1.0;
      } else {
        distf = 1.0 + fabs(gradientx[z][y][x]) * maps[newmapindex + 1] / hx;
      }
    }
  } else {
    if (maps[newmapindex - 1] == -1 && fabs(gradientx[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex - 1] == -1) {
        distf = 1.0;
      } else {
        distf = 1.0 + fabs(gradientx[z][y][x]) * maps[newmapindex - 1] / hx;
      }
    }
  }
  if (gradienty[z][y][x] > 0) {
    if (maps[newmapindex + width] == -1 && fabs(gradienty[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex + width] > -1) {
        distf += fabs(gradienty[z][y][x]) * maps[newmapindex + width] / hy;
      }
    }
  } else {
    if (maps[newmapindex - width] == -1 && fabs(gradienty[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex - width] > -1) {
        distf += fabs(gradienty[z][y][x]) * maps[newmapindex - width] / hy;
      }
    }
  }

  if (gradientz[z][y][x] > 0) {
    if (maps[newmapindex + height * width] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex + height * width] > -1) {
        distf += fabs(gradientz[z][y][x]) * maps[newmapindex + height * width] / hz;
      }
    }
  } else {
    if (maps[newmapindex - height * width] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex - height * width] > -1) {
        distf += fabs(gradientz[z][y][x]) * maps[newmapindex - height * width] / hz;
      }
    }
  }

  distf = distf / (fabs(gradientx[z][y][x]) / hx + fabs(gradienty[z][y][x]) / hy + fabs(gradientz[z][y][x]) / hz);

  return distf;
}

/* Forward-direction counterpart of distanceYezzi_reverse3D; used by
   thickness3DYezzi once the front is past its first, relaxed pass. */
float distanceYezzi3D(float ***gradientx, float ***gradienty, float ***gradientz, int newmapindex, int x, int y, int z, float *maps, int height, int width, float r, float hx, float hy, float hz) {
  float distf;
  if (gradientx[z][y][x] < 0) {
    if (maps[newmapindex + 1] == -1 && fabs(gradientx[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex + 1] == -1) {
        distf = 1.0;
      } else {
        distf = 1.0 + fabs(gradientx[z][y][x]) * maps[newmapindex + 1] / hx;
      }
    }
  } else {
    if (maps[newmapindex - 1] == -1 && fabs(gradientx[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex - 1] == -1) {
        distf = 1.0;
      } else {
        distf = 1.0 + fabs(gradientx[z][y][x]) * maps[newmapindex - 1] / hx;
      }
    }
  }
  if (gradienty[z][y][x] < 0) {
    if (maps[newmapindex + width] == -1 && fabs(gradienty[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex + width] > -1) {
        distf += fabs(gradienty[z][y][x]) * maps[newmapindex + width] / hy;
      }
    }
  } else {
    if (maps[newmapindex - width] == -1 && fabs(gradienty[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex - width] > -1) {
        distf += fabs(gradienty[z][y][x]) * maps[newmapindex - width] / hy;
      }
    }
  }

  if (gradientz[z][y][x] < 0) {
    if (maps[newmapindex + height * width] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex + height * width] > -1) {
        distf += fabs(gradientz[z][y][x]) * maps[newmapindex + height * width] / hz;
      }
    }
  } else {
    if (maps[newmapindex - height * width] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex - height * width] > -1) {
        distf += fabs(gradientz[z][y][x]) * maps[newmapindex - height * width] / hz;
      }
    }
  }

  distf = distf / (fabs(gradientx[z][y][x]) / hx + fabs(gradienty[z][y][x]) / hy + fabs(gradientz[z][y][x]) / hz);
  return distf;
}

/* 3D counterpart of thickness2DYezzi: propagates outward from the voxels
   labeled 0 in `prototypes` (one boundary) across the domain (label 2),
   6-connected, using distanceYezzi3D as the per-step update rule (and the
   relaxed distanceYezzi3D_relax on the first iteration's retry pass, list_aux,
   since few neighbors are populated yet). Runs num_it independent passes;
   `maps` holds the final thickness values on return. */
int thickness3DYezzi(unsigned char *prototypes, int height, int width, int depth, float *maps, float ***gradientx, float ***gradienty, float ***gradientz, int num_it, float hx, float hy, float hz) {
  int i, j, xr, yr, zr, mapindex, newmapindex, d, l, flag, k, count;
  int neighbors[MAX_NEIGHBORS_3D_FACES];
  float distf, r;
  struct index_list list1, list2, list_aux;
  int max_number_in_list = 300000;
  unsigned char *prot_copia;

  prot_copia = (unsigned char *)malloc(sizeof(unsigned char) * height * width * depth);

  for (j = 0; j < height * width * depth; j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list_init(&list1, max_number_in_list);

  list_init(&list2, max_number_in_list);

  list_init(&list_aux, max_number_in_list);
  printf("iteration ");
  for (l = 0; l < num_it; l++) {
    d = 0;
    if (l == 0) {
      r = 0.3;
    } else {
      r = 0.0;
    }
    printf("%02d\b\b", l);
    fflush(0);
    for (i = 0; i < height * width * depth; i++) {
      prot_copia[i] = prototypes[i]; /* to recover the domain after the first iteration */
      if (prot_copia[i] == 0) {
        if (list_push(&list1, i) != 0) {
          printf("Error list1.num_elem >= %d\n", max_number_in_list);
          return 1;
        }
        maps[i] = 0;
      }
    }

    while (list1.num_elem != 0 || list2.num_elem != 0) {
      /* printf("num elem list1: %d\n",list1.num_elem);*/
      while (list1.num_elem != 0) {
        /* Get element from list1 */
        mapindex = list_pop(&list1);
        flag = 0;
        count = neighbors3D_faces(mapindex, height, width, depth, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox3d(newmapindex, height, width);
          yr = maptoy3d(newmapindex, height, width);
          zr = maptoz3d(newmapindex, height, width);
          if (prot_copia[newmapindex] == 2) {
            /* Compute new distance */
            distf = distanceYezzi3D(gradientx, gradienty, gradientz, newmapindex,
                                    xr, yr, zr, maps, height, width, r, hx, hy, hz);

            if (distf > 0 && distf < INF) {
              if (fabs(maps[mapindex] - distf) < 1) {
                maps[newmapindex] = distf;
                flag = 1;
                /* Put new element in list2*/
                if (list_push(&list2, newmapindex) != 0) {
                  printf("Error list2.num_elem >= %d\n", max_number_in_list);
                  return 1;
                }
                prot_copia[newmapindex] = 0;
              }
            }
          }
        }
        if (flag == 0) {
          /* Ponemos en un lista auxiliar el punto que no pudo propagarse*/
          if (list_push(&list_aux, mapindex) != 0) {
            printf("Error list_aux.num_elem >= %d\n", max_number_in_list);
            return 1;
          }
        }
      }

      for (i = 0; i < list_aux.num_elem; i++) {
        /* while (list_aux.num_elem != 0) {    */
        /* Get element from list1 */
        mapindex = list_aux.elem[i];
        count = neighbors3D_faces(mapindex, height, width, depth, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox3d(newmapindex, height, width);
          yr = maptoy3d(newmapindex, height, width);
          zr = maptoz3d(newmapindex, height, width);
          if (prot_copia[newmapindex] == 2) {
            /* Compute new distance */
            if (l == 0) {
              distf = distanceYezzi3D_relax(gradientx, gradienty, gradientz, newmapindex,
                                            xr, yr, zr, maps, height, width, r, hx, hy, hz);
            } else {
              distf = distanceYezzi3D(gradientx, gradienty, gradientz, newmapindex,
                                      xr, yr, zr, maps, height, width, r, hx, hy, hz);
            }
            if (distf > 0 && distf < INF) {
              if (l == 0 || fabs(maps[mapindex] - distf) < 1) {
                maps[newmapindex] = distf;
                /* Put new element in list2*/
                if (list_push(&list2, newmapindex) != 0) {
                  printf("Error list2.num_elem >= %d\n", max_number_in_list);
                  return 1;
                }
                prot_copia[newmapindex] = 0;
              }
            }
          }
        }
      }
      list_clear(&list_aux);
      /* printf("num elem list2: %d\n",list2.num_elem);*/
      d++;
      /* swap(list1.elem,list2.elem); */

      list_swap(&list1, &list2);
    }
  }
  printf("\nmaximum bucket = %d\n", d);

  list_free(&list1);
  list_free(&list2);
  list_free(&list_aux);
  free(prot_copia);
  return 0; /* success */
}

/* Reverse-direction counterpart of thickness3DYezzi: propagates outward from
   the voxels labeled 1 instead of 0, using distanceYezzi_reverse3D (and
   distanceYezzi_reverse3D_relax on the first iteration's retry pass). */
int thickness3DYezzi_reverse(unsigned char *prototypes, int height, int width, int depth, float *maps, float ***gradientx, float ***gradienty, float ***gradientz, int num_it, float hx, float hy, float hz) {
  int i, j, xr, yr, zr, mapindex, newmapindex, d, flag, l, k, count;
  int neighbors[MAX_NEIGHBORS_3D_FACES];
  float distf, r;
  struct index_list list1, list2, list_aux;
  int max_number_in_list = 300000;
  unsigned char *prot_copia;

  prot_copia = (unsigned char *)malloc(sizeof(unsigned char) * height * width * depth);

  for (j = 0; j < height * width * depth; j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list_init(&list1, max_number_in_list);

  list_init(&list2, max_number_in_list);

  list_init(&list_aux, max_number_in_list);
  printf("iteration ");
  for (l = 0; l < num_it; l++) {
    d = 0;
    if (l == 0) {
      r = 0.3;
    } else {
      r = 0.0;
    }
    printf("%02d\b\b", l);
    fflush(0);
    for (i = 0; i < height * width * depth; i++) {
      prot_copia[i] = prototypes[i]; /* to recover the domain after the first iteration */
      if (prot_copia[i] == 1) {
        if (list_push(&list1, i) != 0) {
          printf("Error list1.num_elem >= %d\n", max_number_in_list);
          return 1;
        }
        maps[i] = 0;
      }
    }

    d = 0;
    while (list1.num_elem != 0 || list2.num_elem != 0) {
      /* printf("num elem list1: %d\n",list1.num_elem);*/
      while (list1.num_elem != 0) {
        /* Get element from list1 */
        mapindex = list_pop(&list1);
        flag = 0;
        count = neighbors3D_faces(mapindex, height, width, depth, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox3d(newmapindex, height, width);
          yr = maptoy3d(newmapindex, height, width);
          zr = maptoz3d(newmapindex, height, width);
          if (prot_copia[newmapindex] == 2) {
            /* Compute new distance */
            distf = distanceYezzi_reverse3D(gradientx, gradienty, gradientz, newmapindex,
                                            xr, yr, zr, maps, height, width, r, hx, hy, hz);

            if (distf > 0 && distf < INF) {
              if (fabs(maps[mapindex] - distf) < 1.5) {
                flag = 1;
                maps[newmapindex] = distf;
                /* Put new element in list2*/
                if (list_push(&list2, newmapindex) != 0) {
                  printf("Error list2.num_elem >= %d\n", max_number_in_list);
                  return 1;
                }
                prot_copia[newmapindex] = 1;
              }
            }
          }
        }
        if (flag == 0) {
          /* Ponemos en un lista auxiliar el punto que no pudo propagarse*/
          if (list_push(&list_aux, mapindex) != 0) {
            printf("Error list_aux.num_elem >= %d\n", max_number_in_list);
            return 1;
          }
        }
      }

      for (i = 0; i < list_aux.num_elem; i++) {
        /* while (list_aux.num_elem != 0) { */
        /* Get element from list1 */
        mapindex = list_aux.elem[i];
        count = neighbors3D_faces(mapindex, height, width, depth, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox3d(newmapindex, height, width);
          yr = maptoy3d(newmapindex, height, width);
          zr = maptoz3d(newmapindex, height, width);
          /*if (prototypes[newmapindex] == 2 || prototypes[newmapindex] == 0) {*/
          if (prot_copia[newmapindex] == 2) {
            /* Compute new distance */
            if (l == 0) {
              distf = distanceYezzi_reverse3D_relax(gradientx, gradienty, gradientz, newmapindex,
                                                    xr, yr, zr, maps, height, width, r, hx, hy, hz);
            } else {
              distf = distanceYezzi_reverse3D(gradientx, gradienty, gradientz, newmapindex,
                                              xr, yr, zr, maps, height, width, r, hx, hy, hz);
            }
            if (distf > 0 && distf < INF) {
              if (fabs(maps[mapindex] - distf) < 1.5 || l == 0) {
                maps[newmapindex] = distf;
                /* Put new element in list2*/
                if (list_push(&list2, newmapindex) != 0) {
                  printf("Error list2.num_elem >= %d\n", max_number_in_list);
                  return 1;
                }
                prot_copia[newmapindex] = 1;
              }
            }
          }
        }
      }
      /* printf("num elem list2: %d, list_aux.num_elem %d\n",list2.num_elem,list_aux.num_elem);*/
      list_clear(&list_aux);
      d++;
      /* swap(list1.elem,list2.elem); */

      list_swap(&list1, &list2);
    }
  }
  printf("\nmaximum bucket = %d\n", d);

  list_free(&list1);
  list_free(&list2);
  list_free(&list_aux);
  free(prot_copia);
  return 0; /* success */
}

/* Mean and standard deviation of the thickness map over voxels labeled
   label_cortex that are 6-adjacent to a voxel labeled boundary_l (i.e. the
   band immediately next to one boundary), excluding the outermost 1-voxel
   shell of the volume. *std gets the standard deviation; the mean is returned. */
float compute_mean_thickness(unsigned char *input, float *maps, int label_cortex, int boundary_l, int height, int width, int depth, float *std) {
  int i, j, k, sum = 0, npoints = 0;
  float mean = 0;

  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        if ((i == 0) || (j == 0) || (k == 0) || (i == height - 1) || (j == width - 1) || (k == depth - 1)) {
          /* nothing to do */
        } else if ((input[sum] == label_cortex) &&
                   ((input[sum + 1] == boundary_l) ||
                    (input[sum - 1] == boundary_l) ||

                    (input[sum + height] == boundary_l) ||
                    (input[sum - height] == boundary_l) ||

                    (input[sum + height * width] == boundary_l) ||
                    (input[sum - height * width] == boundary_l))) {
          mean += maps[sum];
          npoints++;
        }
        sum++;
      }
    }
  }

  if (npoints == 0) {
    (*std) = 0;
    return 0;
  }
  mean = mean / (float)npoints;
  sum = 0;
  (*std) = 0;
  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        if ((i == 0) || (j == 0) || (k == 0) || (i == height - 1) || (j == width - 1) || (k == depth - 1)) {
          /* nothing to do */
        } else if ((input[sum] == label_cortex) &&
                   ((input[sum + 1] == boundary_l) ||
                    (input[sum - 1] == boundary_l) ||

                    (input[sum + height] == boundary_l) ||
                    (input[sum - height] == boundary_l) ||

                    (input[sum + height * width] == boundary_l) ||
                    (input[sum - height * width] == boundary_l))) {
          (*std) += (mean - maps[sum]) * (mean - maps[sum]);
        }
        sum++;
      }
    }
  }

  printf("npoints %d\n", npoints);
  (*std) = (npoints > 1) ? sqrt((*std) / (npoints - 1)) : 0;

  return mean;
}

/* Like compute_mean_thickness, but averages over every voxel labeled
   label_cortex in the volume (no boundary-adjacency restriction, no border
   exclusion) rather than just the band next to boundary_l. */
float compute_mean_thickness_volume(unsigned char *input, float *maps, int label_cortex, int height, int width, int depth, float *std) {
  int i, npoints = 0;
  float mean = 0;

  for (i = 0; i < height * width * depth; i++) {
    if (input[i] == label_cortex) {
      mean += maps[i];
      npoints++;
    }
  }
  if (npoints == 0) {
    (*std) = 0;
    return 0;
  }
  mean = mean / (float)npoints;

  (*std) = 0;
  for (i = 0; i < height * width * depth; i++) {
    if (input[i] == label_cortex) {
      (*std) += (mean - maps[i]) * (mean - maps[i]);
    }
  }

  printf("npoints %d\n", npoints);
  (*std) = (npoints > 1) ? sqrt((*std) / (npoints - 1)) : 0;

  return mean;
}
