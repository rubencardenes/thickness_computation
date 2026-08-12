/*(c) Ruben Cardenes Almeida, Valladolid, 8/11/2007 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <poisson2D.h>
#include <laplace2D.h>
#include "png_write.h"

int main(int argc,char* argv[]) {
  unsigned char* input,*min_local;
  float **output;
  int i;
  int iterations = 10;
  float lambda = 0.5,h = 1.0;
  int cortex_label = 2;
  int width, height;
  char *inputfile,*outputfile;
  FILE *fg;

  if (argc != 5 && argc != 4) {
    printf("Usage: poisson2D input.png output.png|output.flt iterations [lambda]\n");
    printf("       (image dimensions are read from the input PNG)\n");
    return 1;
  }

  inputfile  = argv[1];
  outputfile = argv[2];
  iterations = atoi(argv[3]);
  if (argc == 5) {
    lambda = atof(argv[4]);
  }

  /* Read the domain from a PNG; dimensions become height (rows) and width (cols). */
  input = load_png_gray(inputfile, &width, &height);
  if (input == NULL) {
    exit(1);
  }
  printf("Input %s: %d rows x %d cols\n", inputfile, height, width);

  /* Report the distinct label values in the domain. */
  {
    unsigned char present[256];
    print_domain_values(input, height*width, present);
  }

  output = (float**)malloc(sizeof(float*)*height);
  for (i=0;i<height;i++) {
    output[i] = (float*)malloc(sizeof(float)*width);
  }

  printf("Entering in EdgeDetect\n");
  if ( EdgeDetect(input, height, width) == 1 ) {
    printf("Error in EdgeDetect\n");
  }

  printf("Relabeling\n");
  if ( relabel(input,height*width,1,100) != 0) {
    printf("Error in Relabel\n");
  }

  printf("Writing domain domain_anillo_modificado.chr\n");
  fg=fopen("domain_anillo_modificado.chr","w");
  if (fg != NULL) {
    fwrite(input,sizeof(unsigned char),height*width,fg);
    fclose(fg);
  } else {
    fprintf(stderr,"Failed writing domain_anillo_modificado.chr\n");
  }

  printf("Entering in poisson2D\n");
  if ( poisson2D(input, height, width, output, iterations, lambda, 0, h, cortex_label) == 1 ) {
    printf("Error in poisson2D\n");
  }

  min_local = (unsigned char*)malloc(sizeof(unsigned char)*height*width);
  printf("Entering in min local detection\n");
  if ( minimos_locales2D(output, min_local, height, width) == 1 ) {
    printf("Error in minimos_locales2D\n");
  }

  printf("Writing ouput %s:\n",outputfile);
  {
    /* Copy the row-pointer field into a contiguous buffer so it can be written
       as a raw float file or a min..max normalized (optionally colored) PNG. */
    float *buf = (float*)malloc(sizeof(float)*height*width);
    int r,col,idx;
    for (r=0;r<height;r++) {
      for (col=0;col<width;col++) {
        buf[r*width+col] = output[r][col];
      }
    }
    /* For the PNG visualization, show the field only on the solved region
       (input != 0). Background pixels are left non-finite so they render black
       and are excluded from the min..max normalization; the raw .flt output
       keeps the full field. */
    if (png_has_extension(outputfile)) {
      for (idx=0; idx<height*width; idx++) {
        if (input[idx] == 0) buf[idx] = NAN;
      }
    }
    write_float_output(outputfile, buf, width, height, 0, COLOR_GRAY);
    free(buf);
  }

  printf("Writing max_local.png\n");
  {
    /* min_local is 0 everywhere except at local minima, where it holds the
       (truncated) field value; normalize those to 0-255 (0 stays black) so
       the sparse local-minima map is actually visible. */
    float *local_buf = (float*)malloc(sizeof(float)*height*width);
    int idx;
    for (idx=0; idx<height*width; idx++) {
      local_buf[idx] = (float)min_local[idx];
    }
    write_png_from_float("max_local.png", local_buf, width, height, 1, COLOR_GRAY);
    free(local_buf);
  }

  free(min_local);
  free(input);
  free(output);
  printf("OK poisson2D\n");
  return 0;
}
