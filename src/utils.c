#include "utils.h"
#include <stdlib.h>

int maptox(int mapindex, int width) {
  return mapindex % width;
}

int maptoy(int mapindex, int width) {
  return (mapindex - (mapindex % width)) / width;
}

int maptox3d(int mapindex, int height, int width) {
  return (mapindex % (height * width)) % width;
}

int maptoy3d(int mapindex, int height, int width) {
  return ((mapindex % (height * width)) - ((mapindex % (height * width)) % width)) / width;
}

int maptoz3d(int mapindex, int height, int width) {
  return (mapindex - (mapindex % (height * width))) / (height * width);
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

int neighbors3D(int mapindex, int height, int width, int depth, int *neighbors) {
  int plane_size = height * width;
  int plane = mapindex / plane_size;
  int rest = mapindex % plane_size;
  int row = rest / width;
  int col = rest % width;
  int dplane, drow, dcol, nplane, nrow, ncol;
  int count = 0;

  for (dplane = -1; dplane < 2; dplane++) {
    for (dcol = -1; dcol < 2; dcol++) {
      for (drow = -1; drow < 2; drow++) {
        if (dplane == 0 && drow == 0 && dcol == 0) continue;
        nplane = plane + dplane;
        nrow = row + drow;
        ncol = col + dcol;
        if (nplane < 0 || nplane >= depth || nrow < 0 || nrow >= height || ncol < 0 || ncol >= width) continue;
        neighbors[count] = nplane * plane_size + nrow * width + ncol;
        count++;
      }
    }
  }
  return count;
}

int neighbors3D_faces(int mapindex, int height, int width, int depth, int *neighbors) {
  int plane_size = height * width;
  int plane = mapindex / plane_size;
  int rest = mapindex % plane_size;
  int row = rest / width;
  int col = rest % width;
  int dplane, drow, dcol, nplane, nrow, ncol;
  int count = 0;

  for (dplane = -1; dplane < 2; dplane++) {
    for (dcol = -1; dcol < 2; dcol++) {
      for (drow = -1; drow < 2; drow++) {
        if (abs(dplane) + abs(drow) + abs(dcol) != 1) continue;
        nplane = plane + dplane;
        nrow = row + drow;
        ncol = col + dcol;
        if (nplane < 0 || nplane >= depth || nrow < 0 || nrow >= height || ncol < 0 || ncol >= width) continue;
        neighbors[count] = nplane * plane_size + nrow * width + ncol;
        count++;
      }
    }
  }
  return count;
}
