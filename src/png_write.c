/* PNG image I/O for the thickness tools, backed by stb_image / stb_image_write.
 * This is the single translation unit that instantiates the stb implementations,
 * so it must be compiled as C99 or later (see the Makefile). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

#include "png_write.h"

unsigned char *load_png_gray(const char *filename, int *width, int *height) {
  int channels = 1;
  unsigned char *data = stbi_load(filename, width, height, &channels, 1);
  if (data == NULL)
    fprintf(stderr, "Failed reading image %s: %s\n", filename, stbi_failure_reason());
  return data;
}

int png_has_extension(const char *filename) {
  int n = (int) strlen(filename);
  return (n >= 4 && filename[n-4] == '.' &&
          (filename[n-3] == 'p' || filename[n-3] == 'P') &&
          (filename[n-2] == 'n' || filename[n-2] == 'N') &&
          (filename[n-1] == 'g' || filename[n-1] == 'G'));
}

/* --- Color lookup tables (3 rows R,G,B of 256 entries) --- */

static void get_redblue_lut(unsigned char lut[3][256]) {
  int i;
  for (i = 0; i < 256; i++) {
    lut[0][i] = (unsigned char)(255 - i);  /* R */
    lut[1][i] = 0;                          /* G */
    lut[2][i] = (unsigned char) i;          /* B */
  }
}

static void get_random_lut(unsigned char lut[3][256]) {
  int i;
  srand((unsigned) time(NULL));
  for (i = 0; i < 256; i++) {
    lut[0][i] = (unsigned char)(rand() % 256);
    lut[1][i] = (unsigned char)(rand() % 256);
    lut[2][i] = (unsigned char)(rand() % 256);
  }
}

/* Normalize the float map to an 8-bit value per pixel. The thickness algorithms
   can leave non-finite values (e.g. +inf from a divide-by-zero at a boundary
   pixel); those are ignored when finding the range so a single inf/NaN does not
   collapse the image (+inf -> 255, -inf/NaN -> 0). */
static void normalize_to_u8(const float *data, int n, int zero_anchored,
                            unsigned char *out) {
  int i;
  if (zero_anchored) {
    float maxv = 0.0f;
    for (i = 0; i < n; i++)
      if (isfinite(data[i]) && data[i] > maxv) maxv = data[i];
    for (i = 0; i < n; i++) {
      if (!isfinite(data[i])) {
        out[i] = (data[i] > 0.0f) ? 255 : 0;
      } else if (data[i] <= 0.0f || maxv <= 0.0f) {
        out[i] = 0;
      } else {
        int p = (int)(255.0f * data[i] / maxv + 0.5f);
        out[i] = (unsigned char)(p > 255 ? 255 : p);
      }
    }
  } else {
    float mn = 0.0f, mx = 0.0f, range;
    int seen = 0;
    for (i = 0; i < n; i++) {
      if (!isfinite(data[i])) continue;
      if (!seen) { mn = mx = data[i]; seen = 1; }
      else { if (data[i] < mn) mn = data[i]; if (data[i] > mx) mx = data[i]; }
    }
    range = mx - mn;
    for (i = 0; i < n; i++) {
      if (!isfinite(data[i])) {
        out[i] = (data[i] > 0.0f) ? 255 : 0;
      } else if (!seen || range <= 0.0f) {
        out[i] = 0;
      } else {
        int p = (int)(255.0f * (data[i] - mn) / range + 0.5f);
        out[i] = (unsigned char)(p < 0 ? 0 : (p > 255 ? 255 : p));
      }
    }
  }
}

int write_png_from_float(const char *filename, const float *data,
                         int width, int height, int zero_anchored, int color_mode) {
  int n = width * height, i, ok;
  unsigned char *norm;

  if (width <= 0 || height <= 0) return 1;
  norm = (unsigned char *) malloc((size_t) n);
  if (!norm) return 1;
  normalize_to_u8(data, n, zero_anchored, norm);

  if (color_mode == COLOR_GRAY) {
    ok = stbi_write_png(filename, width, height, 1, norm, width);
  } else {
    unsigned char lut[3][256];
    unsigned char *rgb = (unsigned char *) malloc((size_t) n * 3);
    if (!rgb) { free(norm); return 1; }
    if (color_mode == COLOR_RANDOM) get_random_lut(lut);
    else                            get_redblue_lut(lut);  /* default red-blue */
    for (i = 0; i < n; i++) {
      unsigned char v = norm[i];
      if (v == 0) {                 /* background stays black */
        rgb[3*i] = rgb[3*i+1] = rgb[3*i+2] = 0;
      } else {
        rgb[3*i]   = lut[0][v];
        rgb[3*i+1] = lut[1][v];
        rgb[3*i+2] = lut[2][v];
      }
    }
    ok = stbi_write_png(filename, width, height, 3, rgb, width * 3);
    free(rgb);
  }

  free(norm);
  if (!ok) {
    fprintf(stderr, "Failed writing PNG %s\n", filename);
    return 1;
  }
  return 0;
}

int write_float_output(const char *filename, const float *data,
                       int width, int height, int zero_anchored, int color_mode) {
  FILE *fp;

  if (png_has_extension(filename))
    return write_png_from_float(filename, data, width, height, zero_anchored, color_mode);

  fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Failed writing output %s\n", filename);
    return 1;
  }
  fwrite(data, sizeof(float), (size_t) width * height, fp);
  fclose(fp);
  return 0;
}
