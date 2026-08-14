#include "utils.h"
#include <stdlib.h>

int maptox(int mapindex, int width) {
  return mapindex % width;
}

int maptoy(int mapindex, int width) {
  return (mapindex - (mapindex % width)) / width;
}

/* The 3D volumes are laid out i-fastest:

       mapindex = i + j * height + k * height * width

   with i in [0,height), j in [0,width), k in [0,depth) -- the convention
   compute_boundary_cortex3D, edge_detect3D, laplace3D and compute_mean_thickness
   all index with (their +/-1, +/-height and +/-height*width offsets).

   These three used to divide by `width` where the layout uses `height`, which
   is invisible while height == width -- as it is in every phantom under data/ --
   and silently transposes i and j otherwise. */
int maptox3d(int mapindex, int height, int width) {
  return (mapindex % (height * width)) % height;
}

int maptoy3d(int mapindex, int height, int width) {
  return (mapindex % (height * width)) / height;
}

int maptoz3d(int mapindex, int height, int width) {
  return mapindex / (height * width);
}

int mapIndex3D(int r, int c, int z, int nr, int nc, int nz) {
  if (c >= nc) return -1;
  if (c < 0) return -1;
  if (r >= nr) return -1;
  if (r < 0) return -1;
  if (z >= nz) return -1;
  if (z < 0) return -1;
  return c + r * nc + z * nr * nc;
}

void print_timing(FILE *fp, struct timeval start, struct timeval end) {
  double tuend = 1e-06 * (double)end.tv_usec;
  double tustart = 1e-06 * (double)start.tv_usec;
  double tend = end.tv_sec + tuend;
  double tstart = start.tv_sec + tustart;
  fprintf(fp, "Elapsed time: %g\n", (tend - tstart));
}

/* --- index worklist (see utils.h) ---------------------------------------- */

int list_init(struct index_list *list, int max_elem) {
  list->elem = (int *)malloc(sizeof(int) * max_elem);
  list->num_elem = 0;
  list->first = 0;
  list->max_elem = max_elem;
  if (list->elem == NULL) {
    list->max_elem = 0;
    return 1;
  }
  return 0;
}

void list_free(struct index_list *list) {
  free(list->elem);
  list->elem = NULL;
  list->num_elem = 0;
  list->first = 0;
  list->max_elem = 0;
}

void list_clear(struct index_list *list) {
  list->num_elem = 0;
  list->first = 0;
}

int list_push(struct index_list *list, int mapindex) {
  if (list->num_elem >= list->max_elem) {
    return 1;
  }
  list->elem[list->num_elem] = mapindex;
  list->num_elem++;
  return 0;
}

int list_pop(struct index_list *list) {
  if (list->num_elem <= list->first) {
    return -1;
  }
  list->num_elem--;
  return list->elem[list->num_elem];
}

int list_take(struct index_list *list) {
  if (list->first >= list->num_elem) {
    return -1;
  }
  return list->elem[list->first++];
}

void list_swap(struct index_list *a, struct index_list *b) {
  struct index_list tmp = *a;
  *a = *b;
  *b = tmp;
}

/* --- adjacency (see utils.h) --------------------------------------------- */

int neighbors2D(int mapindex, int height, int width, int *neighbors) {
  int row = mapindex / width;
  int col = mapindex % width;
  int drow, dcol, nrow, ncol;
  int count = 0;

  /* column-major, matching the `for (x..) for (y..)` loops this replaced */
  for (dcol = -1; dcol < 2; dcol++) {
    for (drow = -1; drow < 2; drow++) {
      if (drow == 0 && dcol == 0) continue;
      nrow = row + drow;
      ncol = col + dcol;
      if (nrow < 0 || nrow >= height || ncol < 0 || ncol >= width) continue;
      neighbors[count] = nrow * width + ncol;
      count++;
    }
  }
  return count;
}

/* Both 3D versions decompose mapindex with the i-fastest layout documented at
   maptox3d, and step by 1 (i), height (j) and height*width (k).

   The loop nesting -- k outermost, then i, then j -- is not arbitrary: it
   reproduces, offset for offset, the order the previous width-based code
   emitted on a cubic volume, so callers that depend on visit order (the
   thickness front, which may update several neighbours in one call) are
   unaffected where height == width. */
int neighbors3D(int mapindex, int height, int width, int depth, int *neighbors) {
  int plane_size = height * width;
  int k = mapindex / plane_size;
  int rest = mapindex % plane_size;
  int j = rest / height;
  int i = rest % height;
  int di, dj, dk, ni, nj, nk;
  int count = 0;

  for (dk = -1; dk < 2; dk++) {
    for (di = -1; di < 2; di++) {
      for (dj = -1; dj < 2; dj++) {
        if (dk == 0 && di == 0 && dj == 0) continue;
        nk = k + dk;
        nj = j + dj;
        ni = i + di;
        if (nk < 0 || nk >= depth || nj < 0 || nj >= width || ni < 0 || ni >= height) continue;
        neighbors[count] = nk * plane_size + nj * height + ni;
        count++;
      }
    }
  }
  return count;
}

int neighbors3D_faces(int mapindex, int height, int width, int depth, int *neighbors) {
  int plane_size = height * width;
  int k = mapindex / plane_size;
  int rest = mapindex % plane_size;
  int j = rest / height;
  int i = rest % height;
  int di, dj, dk, ni, nj, nk;
  int count = 0;

  for (dk = -1; dk < 2; dk++) {
    for (di = -1; di < 2; di++) {
      for (dj = -1; dj < 2; dj++) {
        if (abs(dk) + abs(di) + abs(dj) != 1) continue;
        nk = k + dk;
        nj = j + dj;
        ni = i + di;
        if (nk < 0 || nk >= depth || nj < 0 || nj >= width || ni < 0 || ni >= height) continue;
        neighbors[count] = nk * plane_size + nj * height + ni;
        count++;
      }
    }
  }
  return count;
}
