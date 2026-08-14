/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include "laplace3D.h"
#include "compute_boundary_cortex3D.h"

int compute_boundary_cortex3D(unsigned short *segmented, int height, int width, int depth, int label_wm, int label_cortex) {
  int i, j, k, sum;

  for (i = 0; i < height * width * depth; i++) {
    if (segmented[i] != label_cortex && segmented[i] != label_wm && segmented[i] != 0) {
      relabel_ushort(segmented, height * width * depth, segmented[i], 0);
    }
  }

  sum = 0;
  for (k = 0; k < depth; k++) {
    for (j = 0; j < width; j++) {
      for (i = 0; i < height; i++) {
        if (k > 0 && j > 0 && i > 0 && k < depth - 1 && j < width - 1 && i < height - 1) {
          if ((segmented[sum] == label_wm) &&
              ((segmented[sum + 1] == label_cortex) ||
               (segmented[sum - 1] == label_cortex) ||

               (segmented[sum + height] == label_cortex) ||
               (segmented[sum - height] == label_cortex) ||

               (segmented[sum + height * width] == label_cortex) ||
               (segmented[sum - height * width] == label_cortex))) {
            segmented[sum] = 1;
          } else if ((segmented[sum] == 0) &&
                     ((segmented[sum + 1] == label_cortex) ||
                      (segmented[sum - 1] == label_cortex) ||

                      (segmented[sum + height] == label_cortex) ||
                      (segmented[sum - height] == label_cortex) ||

                      (segmented[sum + height * width] == label_cortex) ||
                      (segmented[sum - height * width] == label_cortex))) {
            segmented[sum] = 128;
          }
        }
        sum++;
      }
    }
  }

  relabel_ushort(segmented, height * width * depth, 0, 255);
  relabel_ushort(segmented, height * width * depth, label_wm, 255);
  relabel_ushort(segmented, height * width * depth, label_cortex, 2);
  relabel_ushort(segmented, height * width * depth, 128, 0);
  return 0;
}
