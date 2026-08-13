/* Copyright (c) Ruben Cardenes Almeida 08/04/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "thickness2D.h"
#include "laplace2D.h"

#define PI 3.1415927

static int numelembucket[NUM_BUCKETS];
extern int numrechazos;
extern int numasignaciones;
extern int asignacionesraras;
extern int numPrototypes;
int highestIndexClass;
int actualDimension;
int numPrototypesInClass[MAXCLASSNUMBER];
char buffer[2048];
int pdim;

/* Yezzi PDE thickness propagation, reverse direction: estimates the thickness
   at (x,y) from whichever of its x/y neighbors the Laplace gradient points back
   to, weighted by the gradient magnitude. Any already-propagated neighbor is
   used regardless of gradient magnitude.

   `r` is accepted but ignored. The r-tolerance variants that rejected a step
   when the upstream neighbor was unreached and the gradient exceeded `r` were
   dead code and have been removed; the 3D distanceYezzi3D still applies it. So
   in 2D the retry pass in thickness2DYezzi behaves exactly like the main pass,
   whatever `r` the caller computes. */
float distanceYezzi_reverse(float **gradientx, float **gradienty, int newmapindex, int x, int y, float *maps, int width, float r, float hx, float hy) {
  float distf;

  if (gradientx[y][x] > 0) {
    if (maps[newmapindex + 1] == -1) {
      distf = 1.0;
    } else {
      distf = 1.0 + fabs(gradientx[y][x]) * maps[newmapindex + 1] / hx;
    }
  } else {
    if (maps[newmapindex - 1] == -1) {
      distf = 1.0;
    } else {
      distf = 1.0 + fabs(gradientx[y][x]) * maps[newmapindex - 1] / hx;
    }
  }
  if (gradienty[y][x] > 0) {
    if (maps[newmapindex + width] > -1) {
      distf += fabs(gradienty[y][x]) * maps[newmapindex + width] / hy;
    }
  } else {
    if (maps[newmapindex - width] > -1) {
      distf += fabs(gradienty[y][x]) * maps[newmapindex - width] / hy;
    }
  }
  distf = distf / (fabs(gradientx[y][x]) / hx + fabs(gradienty[y][x]) / hy);
  return distf;
}

/* Forward-direction counterpart of distanceYezzi_reverse, following the
   gradient with the opposite sign convention; used together with
   thickness2DYezzi. `r` is accepted but ignored, as above. */
float distanceYezzi(float **gradientx, float **gradienty, int newmapindex, int x, int y, float *maps, int width, float r, float hx, float hy) {
  float distf;
  if (gradientx[y][x] < 0) {
    if (maps[newmapindex + 1] == -1) {
      distf = 1.0;
    } else {
      distf = 1.0 + fabs(gradientx[y][x]) * maps[newmapindex + 1] / hx;
    }
  } else {
    if (maps[newmapindex - 1] == -1) {
      distf = 1.0;
    } else {
      distf = 1.0 + fabs(gradientx[y][x]) * maps[newmapindex - 1] / hx;
    }
  }
  if (gradienty[y][x] < 0) {
    if (maps[newmapindex + width] > -1) {
      distf += fabs(gradienty[y][x]) * maps[newmapindex + width] / hy;
    }
  } else {
    if (maps[newmapindex - width] > -1) {
      distf += fabs(gradienty[y][x]) * maps[newmapindex - width] / hy;
    }
  }
  distf = distf / (fabs(gradientx[y][x]) / hx + fabs(gradienty[y][x]) / hy);
  return distf;
}


/* Computes the 2D thickness map by propagating distf values outward from the
   pixels labeled 0 in `prototypes` (one boundary) across the label_cortex band,
   using distanceYezzi as the per-step update rule. The front is advanced with
   three worklists (list1 = current front, list2 = next front, list_aux = pixels
   that could not propagate to any neighbor this pass and get one more try after
   the front has advanced). Runs num_it independent passes; `maps` holds the
   final thickness values on return.

   The retry pass passes `r` where the main pass passes 0, but distanceYezzi
   ignores `r` in 2D, so the two passes currently apply the same rule. */
int thickness2DYezzi(unsigned char *prototypes, int height, int width, float *maps, float **laplacefield, float **gradientx, float **gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug) {
  int i, j, xr, yr, mapindex, newmapindex, d, l, flag, k, count;
  int neighbors[MAX_NEIGHBORS_2D];
  float distf, r;
  struct index_list list1, list2, list_aux;
  int max_number_in_list = 50000;
  unsigned char *prot_copia;

  prot_copia = (unsigned char *)malloc(sizeof(unsigned char) * height * width);

  for (j = 0; j < height * width; j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list_init(&list1, max_number_in_list);

  list_init(&list2, max_number_in_list);

  list_init(&list_aux, max_number_in_list);

  for (l = 0; l < num_it; l++) {
    d = 0;
    if (l == 0) {
      r = 0.3;
    } else {
      r = 0.04;
    }
    if (debug == 1) {
      printf("iteration %d\n", l);
    }
    for (i = 0; i < height * width; i++) {
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
        count = neighbors2D(mapindex, height, width, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox(newmapindex, width);
          yr = maptoy(newmapindex, width);
          if (prot_copia[newmapindex] == label_cortex) {
            /* Compute new distance */
            distf = distanceYezzi(gradientx, gradienty, newmapindex, xr, yr, maps, width, 0.00, hx, hy);

            if (distf > 0) {
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
        count = neighbors2D(mapindex, height, width, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox(newmapindex, width);
          yr = maptoy(newmapindex, width);
          if (prot_copia[newmapindex] == label_cortex) {
            /* Compute new distance */
            distf = distanceYezzi(gradientx, gradienty, newmapindex, xr, yr, maps, width, r, hx, hy);

            if (distf > 0) {
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
      list_clear(&list_aux);
      /* printf("num elem list2: %d\n",list2.num_elem);*/
      d++;
      /* swap(list1.elem,list2.elem); */

      list_swap(&list1, &list2);
    }
  }
  printf("dmax = %d\n", d);

  list_free(&list1);
  list_free(&list2);
  list_free(&list_aux);
  free(prot_copia);
  return 0; /* success */
}

/* Reverse-direction counterpart of thickness2DYezzi: propagates outward from
   the pixels labeled 1 instead of 0, using distanceYezzi_reverse. Produces the
   thickness map measured from the opposite boundary. */
int thickness2DYezzi_reverse(unsigned char *prototypes, int height, int width, float *maps, float **laplacefield, float **gradientx, float **gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug) {
  int i, j, xr, yr, mapindex, newmapindex, d, flag, l, k, count;
  int neighbors[MAX_NEIGHBORS_2D];
  float distf, r;
  struct index_list list1, list2, list_aux;
  int max_number_in_list = 50000;
  unsigned char *prot_copia;

  prot_copia = (unsigned char *)malloc(sizeof(unsigned char) * height * width);

  for (j = 0; j < height * width; j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list_init(&list1, max_number_in_list);

  list_init(&list2, max_number_in_list);

  list_init(&list_aux, max_number_in_list);

  for (l = 0; l < num_it; l++) {
    d = 0;
    if (l == 0) {
      r = 0.08;
    } else {
      r = 0.00;
    }
    if (debug == 1) {
      printf("iteration %d\n", l);
    }
    for (i = 0; i < height * width; i++) {
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
        count = neighbors2D(mapindex, height, width, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox(newmapindex, width);
          yr = maptoy(newmapindex, width);
          if (prot_copia[newmapindex] == label_cortex) {
            /* Compute new distance */
            distf = distanceYezzi_reverse(gradientx, gradienty, newmapindex, xr, yr, maps, width, 0.000, hx, hy);

            if (distf > 0) {
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
        count = neighbors2D(mapindex, height, width, neighbors);
        for (k = 0; k < count; k++) {
          newmapindex = neighbors[k];
          xr = maptox(newmapindex, width);
          yr = maptoy(newmapindex, width);
          if (prot_copia[newmapindex] == label_cortex) {
            /* Compute new distance */
            distf = distanceYezzi_reverse(gradientx, gradienty, newmapindex, xr, yr, maps, width, r, hx, hy);

            if (distf > 0) {
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

      /* printf("num elem list2: %d, list_aux.num_elem %d\n",list2.num_elem,list_aux.num_elem);*/
      list_clear(&list_aux);
      d++;
      /* swap(list1.elem,list2.elem); */

      list_swap(&list1, &list2);
    }
  }
  printf("dmax = %d\n", d);

  list_free(&list1);
  list_free(&list2);
  list_free(&list_aux);
  free(prot_copia);
  return 0; /* success */
}

/* Quantizes the gradient direction (gradientx,gradienty) at pixel (x,y) into a
   step to one of its 8 neighbors (written to *newx,*newy), carrying forward the
   per-pixel angle_error so the discretization error doesn't accumulate along a
   streamline. Returns the residual angle error for the new pixel. Used by
   thickness2Dgradient to trace paths along the Laplace gradient field. */
float getnewcoordinates(int *newx, int *newy, int x, int y, int width, float gradientx, float gradienty, float *angle_error) {
  double theta, newtheta, dy, dx;
  float new_angle_error = 0;
  theta = acos((double)gradientx);
  if (gradienty < 0) {
    theta = 2 * PI - theta;
  }

  newtheta = theta + angle_error[y * width + x];
  dy = sin(newtheta);
  dx = cos(newtheta);

  if (dy >= 0) {              /* theta -> (0,PI)*/
    if (fabs(dx) <= 0.3826) { /*theta -> (3PI/8,5PI/8) */
      (*newx) = x;
      (*newy) = y - 1;
      new_angle_error = newtheta - PI / 2;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - PI / 2;*/
    }
    if (dx >= 0.9238) { /*theta -> (0,PI/8) */
      (*newx) = x + 1;
      (*newy) = y;
      new_angle_error = newtheta;
      /* angle_error[(*newx)+(*newy)*width] = newtheta;*/
    }
    if (dx <= -0.9238) { /*theta -> (7PI/8,PI) */
      (*newx) = x - 1;
      (*newy) = y;
      new_angle_error = newtheta - PI;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - PI;*/
    }
    if (dx <= 0.9238 && dx >= 0.3826) { /*theta -> (PI/8,3PI/8) */
      (*newx) = x + 1;
      (*newy) = y - 1;
      new_angle_error = newtheta - PI / 4;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - PI / 4; */
    }
    if (dx >= -0.9238 && dx <= -0.3826) { /*theta -> (5PI/8,7PI/8) */
      (*newx) = x - 1;
      (*newy) = y - 1;
      new_angle_error = newtheta - 3 * PI / 4;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - 3*PI / 4; */
    }
  }
  if (dy < 0) {
    if (fabs(dx) <= 0.3826) {
      (*newx) = x;
      (*newy) = y + 1;
      new_angle_error = newtheta - 3 * PI / 2;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - 3*PI / 2; */
    }
    if (dx >= 0.9238) {
      (*newx) = x + 1;
      (*newy) = y;
      new_angle_error = newtheta;
      /* angle_error[(*newx)+(*newy)*width] = newtheta; */
    }
    if (dx <= -0.9238) {
      (*newx) = x - 1;
      (*newy) = y;
      new_angle_error = newtheta - PI;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - PI;*/
    }
    if (dx <= 0.9238 && dx >= 0.3826) {
      (*newx) = x + 1;
      (*newy) = y + 1;
      new_angle_error = newtheta - 7 * PI / 4;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - 7*PI / 4; */
    }
    if (dx >= -0.9238 && dx <= -0.3826) {
      (*newx) = x - 1;
      (*newy) = y + 1;
      new_angle_error = newtheta - 5 * PI / 4;
      /* angle_error[(*newx)+(*newy)*width] = newtheta - 5*PI / 4; */
    }
  }

  return new_angle_error;
}

/* Alternative thickness computation: instead of propagating distance values
   through a wavefront (as thickness2DYezzi does), walks pixel-by-pixel from
   seed points on the label_cortex band along the quantized gradient direction
   (via getnewcoordinates) until it reaches the far boundary, copying the
   corresponding value from a precomputed `input_maps`. */
int thickness2Dgradient(unsigned char *prototypes, int height, int width, float *input_maps,
                        float *maps, float **laplacefield, float **gradientx, float **gradienty) {
  int i, j, x, y, xr, yr, newx, newy, provx, provy, mapindex, newmapindex, provmapindex, aux, d;
  int counter = 0;
  float distf, xf, yf, new_angle_error;
  float *angle_error;
  int numelemaislados = 0;
  struct index_list list1, list2, list3;
  int max_number_in_list = 50000;

  angle_error = (float *)malloc(sizeof(float) * height * width);
  for (j = 0; j < height * width; j++) {
    maps[j] = -1;
    angle_error[j] = 0;
  }

  list_init(&list1, max_number_in_list);

  list_init(&list2, max_number_in_list);

  list_init(&list3, max_number_in_list);

  for (i = 0; i < height * width; i = i + 10) {
    if (prototypes[i] == 2 &&
        (prototypes[i + 1] == 0 || prototypes[i - 1] == 0 ||
         prototypes[i + width] == 0 || prototypes[i - width] == 0)) {
      if (list_push(&list1, i) != 0) {
        printf("Error list1.num_elem >= %d\n", max_number_in_list);
        return 1;
      }
      maps[i] = 0;
    }
  }

  d = 0;
  while (1) {
    while (list1.num_elem != 0) {
      /* printf("num elem list1: %d\n",list1.num_elem);*/
      for (i = 0; i < list1.num_elem; i++) {
        /* Get element from list1 */
        mapindex = list1.elem[i];
        x = maptox(mapindex, width);
        y = maptoy(mapindex, width);
        new_angle_error = getnewcoordinates(&newx, &newy, x, y, width, gradientx[y][x], -gradienty[y][x], angle_error);
        newmapindex = newy * width + newx;

        angle_error[newmapindex] = new_angle_error;
        if (newmapindex == mapindex) {
          printf("newmapindex == mapindex!\n");
        }
        if (prototypes[newmapindex] == 2) {
          /* Compute new dstance */
          /* maps[newmapindex] = maps[mapindex] + distance(x,y,newx,newy);*/
          maps[newmapindex] = input_maps[newmapindex];
          if (maps[newmapindex] > 0) {
            /* Put new element in list2 */
            if (list_push(&list2, newmapindex) != 0) {
              printf("Error list2.num_elem >= %d\n", max_number_in_list);
              return 1;
            }
            /* Mark the new visited site */
            prototypes[newmapindex] = 0;
          }
        }
      }
      /* for (i=0;i<list1.num_elem;i++) {
      mapindex = list1.elem[i];
      for (x=-1;x<2;x++) {  
  for (y=-1;y<2;y++) {
    if (x==0 && y ==0) continue;  
    newmapindex = mapindex + x + y*width;  
    if (prototypes[newmapindex] == 2) {	      
      xr = maptox(mapindex,width);  
      yr = maptoy(mapindex,width);    
      newx = maptox(newmapindex,width);
      newy = maptoy(newmapindex,width);
      new_angle_error = getnewcoordinates(&provx,&provy,newx,newy,width,-gradientx[newy][newx],gradienty[newy][newx],angle_error);
            hay que calcular el angulo de error en (newx,newy)
      provmapindex = provy*width+provx; 
      angle_error[newmapindex]=-new_angle_error;
      if (maps[provmapindex] == -1) {
        maps[newmapindex] = input_maps[newmapindex];
      } else {
        maps[newmapindex] = input_maps[newmapindex];
      }	  
     
      if (maps[newmapindex] != -1) {
        list3.elem[list3.num_elem] = newmapindex; 
        list3.num_elem++;	   
        prototypes[newmapindex] = 0;	  
        numelemaislados++; 
      }
    }
  } 
      } 
    }
    d++;*/
      /* swap(list1.elem,list2.elem); */

      list_swap(&list1, &list2);
      list_clear(&list2);
    }
    /* swap(list1.elem,list3.elem); */
    list_swap(&list1, &list3);
    list_clear(&list3);
    if (list1.num_elem == 0) break;
  }
  printf("Number of collisions %d\n", counter);
  printf("dmax = %d, numelemaislados %d\n", d, numelemaislados);

  free(angle_error);
  list_free(&list1);
  list_free(&list2);
  list_free(&list3);
  return 0; /* success */
}

/* Mean and standard deviation of the thickness over the band (pixels whose
   input label equals label_cortex). Non-finite map values (e.g. +inf left at a
   boundary pixel) are skipped. Returns the mean; *std gets the standard
   deviation and *npoints the number of band pixels averaged. */
float compute_mean_thickness2D(unsigned char *input, float *maps, int label_cortex, int height, int width, int *npoints, float *std) {
  int i, n = 0;
  float mean = 0;

  for (i = 0; i < height * width; i++) {
    if (input[i] == label_cortex && isfinite(maps[i])) {
      mean += maps[i];
      n++;
    }
  }
  if (n == 0) {
    *npoints = 0;
    *std = 0;
    return 0;
  }
  mean = mean / (float)n;

  *std = 0;
  for (i = 0; i < height * width; i++) {
    if (input[i] == label_cortex && isfinite(maps[i])) {
      *std += (mean - maps[i]) * (mean - maps[i]);
    }
  }
  *std = (n > 1) ? sqrt(*std / (n - 1)) : 0;
  *npoints = n;

  return mean;
}
