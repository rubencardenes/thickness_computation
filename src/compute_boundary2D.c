/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "io.h"
#include "laplace2D.h"

int test_boundary(unsigned short *segmented, unsigned short *original, int i, int j, int width, float threshold) {
  int x, y, mapindex;
  int count = 0;
  int sum = 0;
  float media;
  int n_size = 6;
  mapindex = i * width + j;
  for (x = -n_size; x <= n_size; x++) {
    for (y = -n_size; y <= n_size; y++) {
      if (x == 0 && y == 0) continue;
      if (segmented[mapindex + width * y + x] != segmented[mapindex]) {
        count++;
        sum += original[mapindex + width * y + x];
      }
    }
  }

  media = (float)sum / (float)count;
  if (media < threshold) {
    segmented[mapindex] = 1;
  } else {
    segmented[mapindex] = 2;
  }

  return 0;
}

int compute_boundary2D(unsigned short *segmented, unsigned short *original, int height, int width, int label, float threshold) {
  int i, j, istart, jstart, sum;

  istart = -1;
  jstart = -1;
  sum = 0;
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      if (segmented[sum] == label) {
        istart = i;
        jstart = j;
      }
      if (segmented[sum] != 0 && segmented[sum] != label) {
        relabel_ushort(segmented, height * width, segmented[sum], 0);
      }
      sum++;
    }
  }

  if (istart == -1 || jstart == -1) {
    printf("Label %d NOT found \n", label);
    return 0;
  }

  printf("label = %d found istart %d jstart %d \n", label, istart, jstart);
  sum = 0;
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      /* test_boundary scans a window of radius n_size=6, so the margin
   here must be at least 6 to keep that scan in bounds. */
      if ((i < 6) || (j < 6) || (i >= height - 6) || (j >= width - 6)) {
        /* nothing to do */
      } else if ((segmented[sum] == label) &&
                 ((segmented[sum + 1] == 0) ||
                  (segmented[sum - 1] == 0) ||

                  (segmented[sum + width] == 0) ||
                  (segmented[sum - width] == 0))) {
        test_boundary(segmented, original, i, j, width, threshold);
      }
      sum++;
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  unsigned short *segmented;
  unsigned short *original;
  int i, swapbyte, label, max_component, fsize, hsize;
  float threshold;
  int height = 256;
  int width = 256;
  FILE *fp;
  struct timeval startinit;
  char *segmented_file, *original_file;
  unsigned char *aux, *output;

  if (argc != 9) {
    printf("Usage: compute_boundary2D height width segmented_file.ush original_file.mri output2D.chr label threshold swapbyte(0/1) \n");
    return 1;
  }

  height = atoi(argv[1]);
  width = atoi(argv[2]);
  segmented_file = argv[3];
  original_file = argv[4];
  label = atoi(argv[6]);
  threshold = atof(argv[7]);
  swapbyte = atoi(argv[8]);

  gettimeofday(&startinit, NULL);
  /* reserve data */
  segmented = (unsigned short *)malloc(sizeof(unsigned short) * height * width);
  original = (unsigned short *)malloc(sizeof(unsigned short) * height * width);

  /* Read data */
  fp = fopen(segmented_file, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed reading inputfile %s\n", segmented_file);
    exit(1);
  }
  fread(segmented, sizeof(unsigned short), height * width, fp);
  fclose(fp);

  fsize = fileSize(original_file);
  hsize = headerSize(fsize);

  fp = fopen(original_file, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed reading inputfile %s\n", original_file);
    exit(1);
  }
  if (hsize > 0) {
    fseek(fp, hsize, 1);
  }
  fread(original, sizeof(unsigned short), height * width, fp);
  fclose(fp);
  /* If we have to swap the data: */
  if (swapbyte == 1) {
    for (i = 0; i < height * width; i++) {
      aux = (unsigned char *)&segmented[i];
      segmented[i] = ReadGEShort(aux);
      aux = (unsigned char *)&original[i];
      original[i] = ReadGEShort(aux);
    }
  }

  /* Compute */
  compute_boundary2D(segmented, original, height, width, label, threshold);
  sizefilter2D(segmented, height, width, 15, 1, 2);
  sizefilter2D(segmented, height, width, 20, 2, 1);

  /* Keep only the largest label-1 component: measure it, then drop every
     component strictly smaller. sizefilter2D removes components of `max_size`
     pixels or fewer, so the threshold is one below the largest. Both functions
     measure components with collect_component2D, so the two sizes are directly
     comparable -- they must stay that way for this to keep exactly one
     component. */
  max_component = maxcomponent2D(segmented, height, width, 1);
  if (max_component < 0) {
    return 1;
  }
  sizefilter2D(segmented, height, width, max_component - 1, 1, 2);

  relabel_ushort(segmented, height * width, 0, 255);
  relabel_ushort(segmented, height * width, 2, 0);
  relabel_ushort(segmented, height * width, label, 2);

  printf("computing corners\n");
  new_compute_corners(segmented, height, width);

  /* convert data */
  output = (unsigned char *)malloc(height * width * sizeof(unsigned char));
  for (i = 0; i < height * width; i++) {
    output[i] = segmented[i];
  }

  /* Write data */
  fp = fopen(argv[5], "w");
  if (fp == NULL) {
    fprintf(stderr, "Failed writing file %s\n", argv[5]);
    exit(1);
  }
  fwrite(output, sizeof(unsigned char), height * width, fp);

  fclose(fp);

  /* Free data */
  free(original);
  free(segmented);
  free(output);
  printf("compute_boundary2D OK\n");
  return 0;
}
