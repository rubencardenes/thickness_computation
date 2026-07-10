/* PNG image I/O for the thickness tools, backed by stb_image / stb_image_write.
 * Supports reading a grayscale PNG (dimensions taken from the file) and writing
 * a float map as a grayscale or color-mapped PNG, or as a raw float file. */
#ifndef PNG_WRITE_H
#define PNG_WRITE_H

/* Color modes for PNG output (same convention as geodesicDT2d):
     0 = grayscale, 1 = red-blue, 2 = random */
#define COLOR_GRAY   0
#define COLOR_REDBLUE 1
#define COLOR_RANDOM 2

/* Load a PNG as a single-channel (grayscale) image. On success returns a
   malloc'd buffer of (*width)*(*height) bytes (free it with free()) and sets
   *width and *height; returns NULL on failure. */
unsigned char *load_png_gray(const char *filename, int *width, int *height);

/* Returns 1 if filename ends with ".png" (case-insensitive), else 0. */
int png_has_extension(const char *filename);

/* Normalize a contiguous, row-major float map and write it as a PNG.
   width = columns, height = rows.
     zero_anchored != 0: values <= 0 map to 0 (black); the largest positive
                         value maps to 255. Suited to thickness maps.
     zero_anchored == 0: linear min..max maps to 0..255. Suited to fields such
                         as the Laplace solution.
     color_mode: COLOR_GRAY writes 8-bit grayscale; COLOR_REDBLUE / COLOR_RANDOM
                 write an RGB image mapping the normalized value through a color
                 lookup table (value 0 stays black as background).
   Returns 0 on success, non-zero on failure. */
int write_png_from_float(const char *filename, const float *data,
                         int width, int height, int zero_anchored, int color_mode);

/* Write a contiguous, row-major float map either as a PNG (per color_mode, when
   filename ends in ".png") or as a raw little-endian float file (any other
   extension, matching the historical output; color_mode is then ignored).
   width = columns, height = rows. Returns 0 on success, non-zero on failure. */
int write_float_output(const char *filename, const float *data,
                       int width, int height, int zero_anchored, int color_mode);

/* Scan data[0..n-1], record which byte values occur into present[256]
   (each entry set to 0 or 1), print the sorted list of distinct values to
   stdout, and return the number of distinct values. */
int print_domain_values(const unsigned char *data, int n, unsigned char present[256]);

/* Returns 1 if `value` is one of the domain values (present[value] set);
   otherwise prints an error naming the option `optname` and returns 0. */
int require_label(const unsigned char present[256], int value, const char *optname);

/* Same as print_domain_values but for an unsigned short volume (used by the 3D
   tools). Prints every distinct value; present[0..255] is filled for values in
   that range so require_label() can be used with the usual small labels. */
int print_domain_values_ushort(const unsigned short *data, int n, unsigned char present[256]);

#endif
