/* Guards the 3D index helpers against the width/height confusion that hid for
   years behind the fact that every phantom under data/ is a cube.

   The 3D pipeline lays volumes out i-fastest:

       mapindex = i + j * height + k * height * width

   so the neighbour strides are 1 (i), height (j) and height*width (k), and
   maptox3d/maptoy3d/maptoz3d must invert exactly that. Anything that divides
   by `width` instead of `height` agrees with this on a cubic volume and
   silently transposes i and j on every other one -- which is why these checks
   deliberately use height != width != depth. */

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

static int failures = 0;

static void check(int cond, const char *what, int got, int expected) {
  if (!cond) {
    printf("FAIL: %s (got %d, expected %d)\n", what, got, expected);
    failures++;
  }
}

/* Round-trip every voxel of a deliberately non-cubic volume through the
   decomposition helpers. */
static void check_roundtrip(int height, int width, int depth) {
  int i, j, k, m;
  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        m = i + j * height + k * height * width;
        check(maptox3d(m, height, width) == i, "maptox3d", maptox3d(m, height, width), i);
        check(maptoy3d(m, height, width) == j, "maptoy3d", maptoy3d(m, height, width), j);
        check(maptoz3d(m, height, width) == k, "maptoz3d", maptoz3d(m, height, width), k);
        if (failures > 8) return; /* don't spam */
      }
    }
  }
}

/* Every neighbour must be geometrically adjacent under the layout above: the
   6-connected version may only move by one of the three axis strides, and the
   26-connected version must stay within one step on each axis. Neither may
   wrap from i == height-1 into the next j column. */
static void check_neighbours(int height, int width, int depth) {
  int nb[MAX_NEIGHBORS_3D];
  int i, j, k, m, n, count;
  int plane = height * width;

  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        m = i + j * height + k * height * width;

        count = neighbors3D_faces(m, height, width, depth, nb);
        for (n = 0; n < count; n++) {
          int d = nb[n] - m;
          int ok = (d == 1 || d == -1 || d == height || d == -height ||
                    d == plane || d == -plane);
          if (!ok) {
            printf("FAIL: neighbors3D_faces produced offset %d at (i=%d j=%d k=%d); "
                   "allowed are +/-1, +/-%d, +/-%d\n", d, i, j, k, height, plane);
            failures++;
          }
          /* a step along i must not change j or k */
          if ((d == 1 || d == -1) && (maptoy3d(nb[n], height, width) != j ||
                                      maptoz3d(nb[n], height, width) != k)) {
            printf("FAIL: neighbors3D_faces wrapped across a column at (i=%d j=%d k=%d)\n", i, j, k);
            failures++;
          }
          if (failures > 8) return;
        }

        count = neighbors3D(m, height, width, depth, nb);
        for (n = 0; n < count; n++) {
          int di = maptox3d(nb[n], height, width) - i;
          int dj = maptoy3d(nb[n], height, width) - j;
          int dk = maptoz3d(nb[n], height, width) - k;
          if (di < -1 || di > 1 || dj < -1 || dj > 1 || dk < -1 || dk > 1 ||
              (di == 0 && dj == 0 && dk == 0)) {
            printf("FAIL: neighbors3D produced non-adjacent (di=%d dj=%d dk=%d) "
                   "at (i=%d j=%d k=%d)\n", di, dj, dk, i, j, k);
            failures++;
          }
          if (failures > 8) return;
        }
      }
    }
  }
}

/* Interior voxels have all their neighbours; corners have the fewest. */
static void check_counts(int height, int width, int depth) {
  int nb[MAX_NEIGHBORS_3D];
  int mid = (height / 2) + (width / 2) * height + (depth / 2) * height * width;
  int c;

  c = neighbors3D_faces(mid, height, width, depth, nb);
  check(c == 6, "interior voxel should have 6 face neighbours", c, 6);
  c = neighbors3D(mid, height, width, depth, nb);
  check(c == 26, "interior voxel should have 26 neighbours", c, 26);
  c = neighbors3D_faces(0, height, width, depth, nb);
  check(c == 3, "corner voxel should have 3 face neighbours", c, 3);
  c = neighbors3D(0, height, width, depth, nb);
  check(c == 7, "corner voxel should have 7 neighbours", c, 7);
}

int main(void) {
  /* height != width != depth, and all mutually non-dividing, so a width/height
     mix-up cannot accidentally land on a valid index. */
  int height = 5, width = 7, depth = 3;

  printf("checking 3D index helpers on a %dx%dx%d (non-cubic) volume\n", height, width, depth);
  check_roundtrip(height, width, depth);
  check_neighbours(height, width, depth);
  check_counts(height, width, depth);

  /* A cubic volume must keep working too. */
  printf("checking 3D index helpers on a 4x4x4 (cubic) volume\n");
  check_roundtrip(4, 4, 4);
  check_neighbours(4, 4, 4);
  check_counts(4, 4, 4);

  if (failures) {
    printf("check_index3d: %d failure(s)\n", failures);
    return 1;
  }
  printf("check_index3d: OK\n");
  return 0;
}
