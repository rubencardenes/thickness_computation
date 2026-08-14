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

/* ---------------------------------------------------------------------------
   Yezzi PDE thickness propagation, 2D.

   Estimates the thickness at (x,y) from the two upwind neighbours the Laplace
   gradient points back to, weighted by the gradient magnitude. `sign` selects
   the direction: +1 walks the gradient backwards (the "reverse" map, seeded
   from the label-1 boundary), -1 walks it forwards (seeded from label 0). An
   axis' upwind neighbour is the +stride one when sign*g > 0, the -stride one
   otherwise.

   Unlike the 3D rule there is no tolerance and no rejection: whatever upwind
   neighbours have been reached contribute, and the result is always returned.
   The r-tolerance variants that rejected a step when an upstream neighbour was
   unreached were dead code and have been removed, which is also why the retry
   pass in thickness2D_propagate applies exactly the same rule as the main
   pass. */
static float yezzi_step2D(float **gradientx, float **gradienty, int newmapindex,
                          int x, int y, float *maps, int width, int sign,
                          float hx, float hy) {
  float gx = gradientx[y][x];
  float gy = gradienty[y][x];
  float distf = 1.0;
  int up;

  up = (sign * gx > 0) ? newmapindex + 1 : newmapindex - 1;
  if (maps[up] > -1) {
    distf += fabs(gx) * maps[up] / hx;
  }

  up = (sign * gy > 0) ? newmapindex + width : newmapindex - width;
  if (maps[up] > -1) {
    distf += fabs(gy) * maps[up] / hy;
  }

  return distf / (fabs(gx) / hx + fabs(gy) / hy);
}

/* Tries to hand the value at `mapindex` to each of its 8 neighbours that is
   still undecided band (label_cortex), pushing every neighbour it updates onto
   `next` and marking it `mark`. Returns 1 if at least one neighbour was
   updated. Sets *overflow if `next` filled up. */
static int propagate_from2D(int mapindex, unsigned char *prot_copia, float *maps,
                            float **gradientx, float **gradienty, int height, int width,
                            int sign, float hx, float hy, unsigned char label_cortex,
                            unsigned char mark, struct index_list *next, int *overflow) {
  int neighbors[MAX_NEIGHBORS_2D];
  int k, count, newmapindex, xr, yr;
  int moved = 0;
  float distf;

  count = neighbors2D(mapindex, height, width, neighbors);
  for (k = 0; k < count; k++) {
    newmapindex = neighbors[k];
    if (prot_copia[newmapindex] != label_cortex) continue;

    xr = maptox(newmapindex, width);
    yr = maptoy(newmapindex, width);
    distf = yezzi_step2D(gradientx, gradienty, newmapindex, xr, yr, maps, width, sign, hx, hy);
    if (!(distf > 0)) continue;

    maps[newmapindex] = distf;
    moved = 1;
    if (list_push(next, newmapindex) != 0) {
      *overflow = 1;
      return moved;
    }
    prot_copia[newmapindex] = mark;
  }
  return moved;
}

/* Computes the 2D thickness map by propagating outward from every pixel labeled
   `seed_label` (one boundary) across the label_cortex band, 8-connected,
   running num_it independent passes. `maps` holds the final values on return.

   Three worklists advance the front: list1 is the current front, list2 the next
   one, and list_aux collects the pixels that could not reach any neighbour this
   round and get one more try once the front has moved. */
static int thickness2D_propagate(unsigned char *prototypes, int height, int width, float *maps,
                                 float **gradientx, float **gradienty, int num_it,
                                 float hx, float hy, unsigned char label_cortex, int debug,
                                 unsigned char seed_label, int sign) {
  int npixels = height * width;
  int max_number_in_list = 50000;
  int i, j, l, d = 0, mapindex;
  int overflow = 0;
  struct index_list list1, list2, list_aux;
  unsigned char *prot_copia;

  prot_copia = (unsigned char *)malloc(sizeof(unsigned char) * npixels);
  if (prot_copia == NULL) {
    fprintf(stderr, "thickness2D: out of memory\n");
    return 1;
  }
  for (j = 0; j < npixels; j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  if (list_init(&list1, max_number_in_list) != 0 ||
      list_init(&list2, max_number_in_list) != 0 ||
      list_init(&list_aux, max_number_in_list) != 0) {
    fprintf(stderr, "thickness2D: out of memory\n");
    list_free(&list1);
    list_free(&list2);
    list_free(&list_aux);
    free(prot_copia);
    return 1;
  }

  for (l = 0; l < num_it && !overflow; l++) {
    d = 0;
    if (debug == 1) {
      printf("iteration %d\n", l);
    }

    /* Reseed: recover the domain, since the previous pass overwrote it. */
    for (i = 0; i < npixels; i++) {
      prot_copia[i] = prototypes[i];
      if (prot_copia[i] == seed_label) {
        if (list_push(&list1, i) != 0) {
          overflow = 1;
          break;
        }
        maps[i] = 0;
      }
    }

    while (!overflow && (list1.num_elem != 0 || list2.num_elem != 0)) {
      while (list1.num_elem != 0) {
        mapindex = list_pop(&list1);
        if (!propagate_from2D(mapindex, prot_copia, maps, gradientx, gradienty,
                              height, width, sign, hx, hy, label_cortex,
                              seed_label, &list2, &overflow)) {
          /* Could not reach anyone; give it one more try after the front moves. */
          if (list_push(&list_aux, mapindex) != 0) {
            overflow = 1;
            break;
          }
        }
        if (overflow) break;
      }

      for (i = 0; i < list_aux.num_elem && !overflow; i++) {
        propagate_from2D(list_aux.elem[i], prot_copia, maps, gradientx, gradienty,
                         height, width, sign, hx, hy, label_cortex,
                         seed_label, &list2, &overflow);
      }
      list_clear(&list_aux);
      d++;
      list_swap(&list1, &list2);
    }
  }

  if (overflow) {
    printf("Error: worklist exceeded %d entries\n", max_number_in_list);
  } else {
    printf("dmax = %d\n", d);
  }

  list_free(&list1);
  list_free(&list2);
  list_free(&list_aux);
  free(prot_copia);
  return overflow ? 1 : 0;
}

/* Propagates outward from the pixels labeled 0 (one boundary), following the
   Laplace gradient forwards. */
int thickness2DYezzi(unsigned char *prototypes, int height, int width, float *maps, float **gradientx, float **gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug) {
  return thickness2D_propagate(prototypes, height, width, maps, gradientx, gradienty,
                               num_it, hx, hy, label_cortex, debug, 0, -1);
}

/* Reverse counterpart: propagates from the pixels labeled 1, walking the
   gradient backwards, giving the thickness measured from the other boundary. */
int thickness2DYezzi_reverse(unsigned char *prototypes, int height, int width, float *maps, float **gradientx, float **gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug) {
  return thickness2D_propagate(prototypes, height, width, maps, gradientx, gradienty,
                               num_it, hx, hy, label_cortex, debug, 1, 1);
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
                        float *maps, float **gradientx, float **gradienty) {
  int i, j, x, y, newx, newy, mapindex, newmapindex, d;
  int counter = 0;
  float new_angle_error;
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
   input label equals label_cortex). Returns the mean; *std gets the standard
   deviation and *npoints the number of band pixels averaged.

   Pixels the front never reached are skipped and reported, not averaged in:
   they still hold the -1 sentinel when the caller has not summed the two
   directions, and averaging that in drags the mean below zero-thickness. The
   same test also catches the +inf a boundary pixel can be left at when the
   gradient there is exactly zero. */
float compute_mean_thickness2D(unsigned char *input, float *maps, int label_cortex, int height, int width, int *npoints, float *std) {
  int i, n = 0, nskipped = 0;
  float mean = 0;

  for (i = 0; i < height * width; i++) {
    if (input[i] == label_cortex) {
      if (isfinite(maps[i]) && maps[i] > 0) {
        mean += maps[i];
        n++;
      } else {
        nskipped++;
      }
    }
  }
  if (nskipped > 0) {
    printf("WARNING: %d of %d band pixels were never reached by the propagation "
           "and are excluded; %.1f%% coverage\n",
           nskipped, n + nskipped, 100.0 * n / (n + nskipped));
  }
  if (n == 0) {
    *npoints = 0;
    *std = 0;
    return 0;
  }
  mean = mean / (float)n;

  *std = 0;
  for (i = 0; i < height * width; i++) {
    if (input[i] == label_cortex && isfinite(maps[i]) && maps[i] > 0) {
      *std += (mean - maps[i]) * (mean - maps[i]);
    }
  }
  *std = (n > 1) ? sqrt(*std / (n - 1)) : 0;
  *npoints = n;

  return mean;
}
