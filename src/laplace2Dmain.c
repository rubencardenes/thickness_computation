/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <getopt.h>
#include "png_write.h"

int main(int argc,char* argv[]) {
  unsigned char* input;
  float **output;
  int i;
  int max1 = 256;
  int max2 = 256;
  int iterations = 10;
  float lambda = 0.5;
  int color_mode = COLOR_GRAY;
  int width, height, c, nargs;
  char *inputfile,*outputfile;
  FILE *fg;

  while ((c = getopt(argc, argv, "c:")) != -1) {
    switch (c) {
    case 'c': color_mode = atoi(optarg); break;
    default: break;
    }
  }

  nargs = argc - optind;
  if (nargs != 3 && nargs != 4) {
    printf("Usage: laplace2D [-c color] input.png output.png iterations [lambda]\n");
    printf("       (image dimensions are read from the input PNG)\n");
    printf("       -c color output (0: gray, 1: red-blue, 2: random)\n");
    return 1;
  }

  inputfile  = argv[optind + 0];
  outputfile = argv[optind + 1];
  iterations = atoi(argv[optind + 2]);
  if (nargs == 4) {
    lambda = atof(argv[optind + 3]);
  }

  /* Read the domain from a PNG; dimensions become max1 (rows) and max2 (cols). */
  input = load_png_gray(inputfile, &width, &height);
  if (input == NULL) {
    exit(1);
  }
  max2 = width;
  max1 = height;
  printf("Input %s: %d rows x %d cols\n", inputfile, max1, max2);

  output = (float**)malloc(sizeof(float*)*max1);
  for (i=0;i<max1;i++) {
    output[i] = (float*)malloc(sizeof(float)*max2);
  }
  
  printf("Entering in EdgeDetect\n");
  if ( EdgeDetect(input, max1, max2) == 1 ) {
    printf("Error in EdgeDetect\n");
  } 
  
  printf("Relabeling\n");
  if ( relabel(input,max1*max2,0,2) != 0) {
    printf("Error in Relabel\n");
  }
  
  printf("Entering in RelabelBoundary\n");
  if ( RelabelBoundary(input, max1, max2) == 1 ) {
    printf("Error in RelabelBoundary\n");
  } 
  
  printf("Writing domain domain_anillo_modificado.chr\n");
  fg=fopen("domain_anillo_modificado.chr","w");
  if (fg != NULL) {
    fwrite(input,sizeof(unsigned char),max1*max2,fg);
    fclose(fg);
  } else {
    fprintf(stderr,"Failed writing domain_anillo_modificado.chr\n");
  }

  printf("Entering in laplacian2D\n");
  if ( laplace2D(input, max1, max2, output, iterations, lambda) == 1 ) {
    printf("Error in thickness2D\n");
  }

  printf("Writing ouput %s:\n",outputfile);
  {
    /* Copy the row-pointer field into a contiguous buffer so it can be written
       as a raw float file or a min..max normalized (optionally colored) PNG. */
    float *buf = (float*)malloc(sizeof(float)*max1*max2);
    int r,col;
    for (r=0;r<max1;r++) {
      for (col=0;col<max2;col++) {
        buf[r*max2+col] = output[r][col];
      }
    }
    write_float_output(outputfile, buf, max2, max1, 0, color_mode);
    free(buf);
  }

  free(input);
  free(output);
  printf("OK laplace2D\n");
  return 0;
}

