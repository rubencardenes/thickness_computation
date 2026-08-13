/* Verifies the adjacency fix: neighbours must be geometrically adjacent, and
   two components separated only by the image edge must not be merged. */
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "laplace2D.h"

static int failures = 0;

static void check(int cond, const char *what) {
  printf("%-62s %s\n", what, cond ? "ok" : "FAIL");
  if (!cond) failures++;
}

int main(void) {
  int neighbors[MAX_NEIGHBORS_2D];
  int height = 6, width = 5;
  int i, k, count, row, col, nrow, ncol, bad, wraps;
  unsigned short *img;

  /* 1. every neighbour of every pixel is geometrically adjacent */
  bad = 0;
  wraps = 0;
  for (i = 0; i < height * width; i++) {
    row = i / width;
    col = i % width;
    count = neighbors2D(i, height, width, neighbors);
    for (k = 0; k < count; k++) {
      nrow = neighbors[k] / width;
      ncol = neighbors[k] % width;
      if (abs(nrow - row) > 1 || abs(ncol - col) > 1) bad++;
      /* the old bug: column 0 picking up column width-1 */
      if (col == 0 && ncol == width - 1) wraps++;
      if (col == width - 1 && ncol == 0) wraps++;
    }
  }
  check(bad == 0, "every neighbour is within one row/column");
  check(wraps == 0, "no neighbour wraps around a row edge");

  /* corner has 3 neighbours, edge 5, interior 8 */
  check(neighbors2D(0, height, width, neighbors) == 3, "corner pixel has 3 neighbours");
  check(neighbors2D(width, height, width, neighbors) == 5, "left-edge pixel has 5 neighbours");
  check(neighbors2D(width + 1, height, width, neighbors) == 8, "interior pixel has 8 neighbours");

  /* 2. two blobs touching opposite edges on adjacent rows must stay separate.
        Row 1 column 0, and row 0 column width-1: under the old offset-based
        rule these were "adjacent" and merged into one component of 2. */
  img = (unsigned short *)calloc(height * width, sizeof(unsigned short));
  img[0 * width + (width - 1)] = 1; /* row 0, last column */
  img[1 * width + 0] = 1;           /* row 1, first column */
  check(maxcomponent2D(img, height, width, 1) == 1,
        "edge-separated pixels are two components of size 1, not one of 2");

  /* 3. a genuinely connected diagonal pair is still one component of 2 */
  free(img);
  img = (unsigned short *)calloc(height * width, sizeof(unsigned short));
  img[1 * width + 1] = 1;
  img[2 * width + 2] = 1;
  check(maxcomponent2D(img, height, width, 1) == 2,
        "diagonally touching pixels are still one component of 2");

  /* 4. sizefilter2D must agree with maxcomponent2D: keep only the largest */
  free(img);
  img = (unsigned short *)calloc(height * width, sizeof(unsigned short));
  img[0] = 1;                 /* component of 1 */
  img[3 * width + 1] = 1;     /* component of 3 */
  img[3 * width + 2] = 1;
  img[4 * width + 2] = 1;
  {
    int max_size = maxcomponent2D(img, height, width, 1);
    int survivors = 0;
    check(max_size == 3, "maxcomponent2D finds the size-3 component");
    sizefilter2D(img, height, width, max_size - 1, 1, 2);
    for (i = 0; i < height * width; i++)
      if (img[i] == 1) survivors++;
    check(survivors == 3, "sizefilter2D(max-1) leaves exactly the largest component");
  }
  free(img);

  printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "all checks passed", failures);
  return failures != 0;
}
