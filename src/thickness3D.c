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

/* ---------------------------------------------------------------------------
   Yezzi PDE thickness propagation, 3D.

   One update rule serves the four variants that used to be written out
   separately (forward/reverse x relaxed/strict). They differ in two things:

     `sign`   +1 walks the gradient backwards -- the "reverse" map, seeded from
              the label-1 boundary; -1 walks it forwards, seeded from label 0.
              An axis' upwind neighbour is the +stride one when sign*g > 0 and
              the -stride one when sign*g < 0.

     `strict` 0 selects the relaxed rule used on the first pass, when few
              neighbours are populated yet: take whatever upwind neighbours have
              been reached, and reject (-1) only if none of the three axes had
              one. 1 selects the rule used afterwards: reject as soon as an
              axis' upwind neighbour is unreached and its gradient exceeds the
              tolerance `r`, which is what stops the front from running ahead of
              its own data.

   An axis whose gradient component is exactly 0 is skipped: no transport along
   it, so its neighbours carry no information. Either way it would contribute
   fabs(0)*value == 0 to the numerator and 0 to the denominator; skipping
   matters only under the relaxed rule, where counting it would accept a voxel
   on the strength of an axis that said nothing. Under the strict rule it cannot
   trigger a rejection either, since fabs(0) > r is false for the non-negative r
   the drivers pass, so skipping is a no-op there. */
static float yezzi_step3D(float ***gradientx, float ***gradienty, float ***gradientz,
                          int newmapindex, int i, int j, int k, float *maps,
                          int height, int width, int sign, int strict, float r,
                          float hx, float hy, float hz) {
  float grad[3], h[3];
  int stride[3];
  float distf = 1.0;
  int reached = 0;
  int a, up;

  /* Which gradient array holds which derivative is not what the names suggest.
     The field and the gradients are allocated [depth][height][width], i.e.
     indexed [k][i][j]; iGradY3D differentiates the middle index (i) and
     iGradX3D the last one (j). So:

         axis i -> gradienty, spacing hy, stride 1
         axis j -> gradientx, spacing hx, stride height
         axis k -> gradientz, spacing hz, stride height*width

     which also matches how laplace3D_voxelsize pairs hx/hy/hz with the axes.
     This used to read gradientx[k][j][i] and pair it with stride 1, so the
     j-derivative drove motion along i and the value was sampled at the i<->j
     mirrored voxel. Both are invisible when the volume and the geometry are
     symmetric in i and j, which is why the sphere phantom never showed it. */
  grad[0] = gradienty[k][i][j];
  grad[1] = gradientx[k][i][j];
  grad[2] = gradientz[k][i][j];
  h[0] = hy;
  h[1] = hx;
  h[2] = hz;
  stride[0] = 1;
  stride[1] = height;
  stride[2] = height * width;

  for (a = 0; a < 3; a++) {
    if (sign * grad[a] > 0) {
      up = newmapindex + stride[a];
    } else if (sign * grad[a] < 0) {
      up = newmapindex - stride[a];
    } else {
      continue;
    }
    if (maps[up] > -1) {
      distf += fabs(grad[a]) * maps[up] / h[a];
      reached = 1;
    } else if (strict && fabs(grad[a]) > r) {
      return -1;
    }
  }

  if (!strict && !reached) {
    return -1;
  }
  return distf / (fabs(grad[0]) / h[0] + fabs(grad[1]) / h[1] + fabs(grad[2]) / h[2]);
}

/* Tries to hand the value at `mapindex` to each of its 6 face neighbours that
   is still undecided domain (label 2), pushing every neighbour it updates onto
   `next` and marking it `mark`. Returns 1 if at least one neighbour was
   updated, which is what tells the caller this voxel does not need retrying.

   `relaxed` selects the first-pass update rule (see yezzi_step3D) and also
   drops the `tol` continuity check, which nothing can satisfy before the map
   has been seeded. Sets *overflow if `next` filled up. */
static int propagate_from3D(int mapindex, unsigned char *prot_copia, float *maps,
                            float ***gradientx, float ***gradienty, float ***gradientz,
                            int height, int width, int depth, int sign, int relaxed,
                            float r, float tol, float hx, float hy, float hz,
                            unsigned char mark, struct index_list *next, int *overflow) {
  int neighbors[MAX_NEIGHBORS_3D_FACES];
  int n, count, newmapindex, ir, jr, kr;
  int moved = 0;
  float distf;

  count = neighbors3D_faces(mapindex, height, width, depth, neighbors);
  for (n = 0; n < count; n++) {
    newmapindex = neighbors[n];
    if (prot_copia[newmapindex] != 2) continue;

    ir = maptox3d(newmapindex, height, width);
    jr = maptoy3d(newmapindex, height, width);
    kr = maptoz3d(newmapindex, height, width);
    distf = yezzi_step3D(gradientx, gradienty, gradientz, newmapindex, ir, jr, kr,
                         maps, height, width, sign, !relaxed, r, hx, hy, hz);

    if (!(distf > 0 && distf < INF)) continue;
    if (!relaxed && !(fabs(maps[mapindex] - distf) < tol)) continue;

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

/* Propagates the thickness front outward from every voxel labeled `seed_label`
   across the domain (label 2), 6-connected, running num_it independent passes.
   `maps` holds the final thickness values on return.

   Three worklists advance the front: list1 is the current front, list2 the next
   one, and list_aux collects the voxels that could not reach any neighbour this
   round and get one more try once the front has moved. That retry uses the
   relaxed rule on the first iteration only.

   `sign` selects the direction (see yezzi_step3D). `tol` bounds how far a newly
   computed thickness may sit from the value it propagated from; the forward and
   reverse directions were tuned with different bounds. */
static int thickness3D_propagate(unsigned char *prototypes, int height, int width, int depth,
                                 float *maps, float ***gradientx, float ***gradienty,
                                 float ***gradientz, int num_it, float hx, float hy, float hz,
                                 unsigned char seed_label, int sign, float tol) {
  int nvoxels = height * width * depth;
  int max_number_in_list = 300000;
  int i, j, l, d = 0, mapindex;
  int overflow = 0;
  float r;
  struct index_list list1, list2, list_aux;
  unsigned char *prot_copia;

  prot_copia = (unsigned char *)malloc(sizeof(unsigned char) * nvoxels);
  if (prot_copia == NULL) {
    fprintf(stderr, "thickness3D: out of memory\n");
    return 1;
  }
  for (j = 0; j < nvoxels; j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  if (list_init(&list1, max_number_in_list) != 0 ||
      list_init(&list2, max_number_in_list) != 0 ||
      list_init(&list_aux, max_number_in_list) != 0) {
    fprintf(stderr, "thickness3D: out of memory\n");
    list_free(&list1);
    list_free(&list2);
    list_free(&list_aux);
    free(prot_copia);
    return 1;
  }

  printf("iteration ");
  for (l = 0; l < num_it && !overflow; l++) {
    d = 0;
    r = (l == 0) ? 0.3 : 0.0;
    printf("%02d\b\b", l);
    fflush(0);

    /* Reseed: recover the domain, since the previous pass overwrote it. */
    for (i = 0; i < nvoxels; i++) {
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
        if (!propagate_from3D(mapindex, prot_copia, maps, gradientx, gradienty, gradientz,
                              height, width, depth, sign, 0, r, tol, hx, hy, hz,
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
        propagate_from3D(list_aux.elem[i], prot_copia, maps, gradientx, gradienty, gradientz,
                         height, width, depth, sign, (l == 0), r, tol, hx, hy, hz,
                         seed_label, &list2, &overflow);
      }
      list_clear(&list_aux);
      d++;
      list_swap(&list1, &list2);
    }
  }

  if (overflow) {
    printf("\nError: worklist exceeded %d entries\n", max_number_in_list);
  } else {
    printf("\nmaximum bucket = %d\n", d);
  }

  list_free(&list1);
  list_free(&list2);
  list_free(&list_aux);
  free(prot_copia);
  return overflow ? 1 : 0;
}

/* Propagates outward from the voxels labeled 0 (one boundary), following the
   Laplace gradient forwards. */
int thickness3DYezzi(unsigned char *prototypes, int height, int width, int depth, float *maps, float ***gradientx, float ***gradienty, float ***gradientz, int num_it, float hx, float hy, float hz) {
  return thickness3D_propagate(prototypes, height, width, depth, maps,
                               gradientx, gradienty, gradientz, num_it, hx, hy, hz,
                               0, -1, 1.0);
}

/* Reverse counterpart: propagates from the voxels labeled 1, walking the
   gradient backwards, and was tuned with a looser continuity bound. */
int thickness3DYezzi_reverse(unsigned char *prototypes, int height, int width, int depth, float *maps, float ***gradientx, float ***gradienty, float ***gradientz, int num_it, float hx, float hy, float hz) {
  return thickness3D_propagate(prototypes, height, width, depth, maps,
                               gradientx, gradienty, gradientz, num_it, hx, hy, hz,
                               1, 1, 1.5);
}

/* True for voxels in the band that touch a `boundary_l` voxel across one of
   their 6 faces, excluding the outermost 1-voxel shell (whose neighbour offsets
   would leave the volume). Both passes of compute_mean_thickness must select
   exactly the same voxels, so they share this one predicate. */
static int is_band_face_adjacent(const unsigned char *input, int sum, int i, int j, int k,
                                 int height, int width, int depth,
                                 int label_cortex, int boundary_l) {
  if ((i == 0) || (j == 0) || (k == 0) || (i == height - 1) || (j == width - 1) || (k == depth - 1)) {
    return 0;
  }
  if (input[sum] != label_cortex) {
    return 0;
  }
  return (input[sum + 1] == boundary_l) ||
         (input[sum - 1] == boundary_l) ||
         (input[sum + height] == boundary_l) ||
         (input[sum - height] == boundary_l) ||
         (input[sum + height * width] == boundary_l) ||
         (input[sum - height * width] == boundary_l);
}

/* A voxel the front never reached carries no thickness and must not be averaged
   as one. The propagation leaves such voxels at the -1 sentinel, but callers
   relabel -1 to 0 before writing the volume out, so by the time the statistics
   run the sentinel is already gone -- hence the test is "not a positive, finite
   value" rather than "not -1". A real thickness is at least one propagation
   step, so no legitimate measurement is excluded by this. */
static int has_thickness(float v) {
  return isfinite(v) && v > 0;
}

/* Print how much of the measured region actually carried a thickness. A large
   skipped count means the front stalled, which is a result worth seeing rather
   than one to average away. */
static void report_coverage(int npoints, int nskipped) {
  int total = npoints + nskipped;
  printf("npoints %d", npoints);
  if (nskipped > 0) {
    printf(" (WARNING: %d of %d measured voxels were never reached by the "
           "propagation and are excluded; %.1f%% coverage)",
           nskipped, total, total ? 100.0 * npoints / total : 0.0);
  }
  printf("\n");
}

/* Mean and standard deviation of the thickness map over voxels labeled
   label_cortex that are 6-adjacent to a voxel labeled boundary_l (i.e. the
   band immediately next to one boundary), excluding the outermost 1-voxel
   shell of the volume. *std gets the standard deviation; the mean is returned.

   Voxels the front never reached are skipped and reported rather than counted
   as zero thickness: averaging them in used to turn "the propagation failed
   over half this surface" into "the tissue is thin", silently. On the ellipsoid
   phantom that was the difference between a reported 4.03 and an actual 8.94. */
float compute_mean_thickness(unsigned char *input, float *maps, int label_cortex, int boundary_l, int height, int width, int depth, float *std) {
  int i, j, k, sum = 0, npoints = 0, nskipped = 0;
  float mean = 0;

  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        if (is_band_face_adjacent(input, sum, i, j, k, height, width, depth, label_cortex, boundary_l)) {
          if (has_thickness(maps[sum])) {
            mean += maps[sum];
            npoints++;
          } else {
            nskipped++;
          }
        }
        sum++;
      }
    }
  }

  report_coverage(npoints, nskipped);
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
        if (is_band_face_adjacent(input, sum, i, j, k, height, width, depth, label_cortex, boundary_l) &&
            has_thickness(maps[sum])) {
          (*std) += (mean - maps[sum]) * (mean - maps[sum]);
        }
        sum++;
      }
    }
  }

  (*std) = (npoints > 1) ? sqrt((*std) / (npoints - 1)) : 0;

  return mean;
}

/* Like compute_mean_thickness, but averages over every voxel labeled
   label_cortex in the volume (no boundary-adjacency restriction, no border
   exclusion) rather than just the band next to boundary_l. Skips and reports
   unreached voxels for the same reason. */
float compute_mean_thickness_volume(unsigned char *input, float *maps, int label_cortex, int height, int width, int depth, float *std) {
  int i, npoints = 0, nskipped = 0;
  float mean = 0;

  for (i = 0; i < height * width * depth; i++) {
    if (input[i] == label_cortex) {
      if (has_thickness(maps[i])) {
        mean += maps[i];
        npoints++;
      } else {
        nskipped++;
      }
    }
  }
  report_coverage(npoints, nskipped);
  if (npoints == 0) {
    (*std) = 0;
    return 0;
  }
  mean = mean / (float)npoints;

  (*std) = 0;
  for (i = 0; i < height * width * depth; i++) {
    if (input[i] == label_cortex && has_thickness(maps[i])) {
      (*std) += (mean - maps[i]) * (mean - maps[i]);
    }
  }

  (*std) = (npoints > 1) ? sqrt((*std) / (npoints - 1)) : 0;

  return mean;
}
