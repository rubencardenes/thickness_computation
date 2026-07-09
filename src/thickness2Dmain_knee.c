/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include "thickness2D.h"
#include "laplace2D.h"
#include "png_write.h"

int numasignaciones = 0;

int main(int argc, char* argv[]) {
  unsigned char *input;
  unsigned char label = 2;
  float *maps;
  float **laplacefield,**gradientx,**gradienty;
  float lambda = 0.5, hx = 1,hy = 1;
  int i,col,row,iterations;
  int max1 = 256;
  int max2 = 256;
  int max3 = 1;
  int reverse = 0;
  int num_it;
  int debug = 1;
  int color_mode = COLOR_GRAY;
  int width, height, c, nargs;
  FILE *fp,*fg;
  struct timeval startinit;
  struct timeval endinit;
  struct timeval endtotal;
  char *inputfile,*outputfile,*laplacefile;

  while ((c = getopt(argc, argv, "dc:")) != -1) {
    switch (c) {
    case 'd': debug = 1; break;
    case 'c': color_mode = atoi(optarg); break;
    default: break;
    }
  }

  nargs = argc - optind;
  if (nargs < 6 || nargs > 8) {
    printf("Usage: thickness2D_knee [options] domain2D.png output2D.png iterations_laplace iterations hx hy [reverse] [lambda]\n");
    printf("       (image dimensions are read from the input PNG)\n");
    printf("       -c color output (0: gray, 1: red-blue, 2: random)\n");
    return 1;
  }

  inputfile  = argv[optind + 0];
  outputfile = argv[optind + 1];
  iterations = atoi(argv[optind + 2]);
  num_it     = atoi(argv[optind + 3]);
  hx = atof(argv[optind + 4]);
  hy = atof(argv[optind + 5]);
  if (nargs >= 7) reverse = atoi(argv[optind + 6]);
  if (nargs == 8) lambda = atof(argv[optind + 7]);

  gettimeofday(&startinit,NULL);

  /* Read the domain from a PNG; dimensions become max1 (rows) and max2 (cols). */
  input = load_png_gray(inputfile, &width, &height);
  if (input == NULL) {
    exit(1);
  }
  max2 = width;
  max1 = height;
  printf("Input %s: %d rows x %d cols\n", inputfile, max1, max2);

  laplacefield  = (float**)malloc(sizeof(float*)*max1);
  for (i=0;i<max1;i++) {
    laplacefield[i] = (float*)malloc(sizeof(float)*max2);
  }

  printf("Entering in laplacian2D\n");
  if ( laplace2D(input, max1, max2, laplacefield, iterations, lambda, reverse) == 1 ) {
    printf("Error in thickness2D\n");
  }

  if (debug == 1) {
    printf("Writing ouput input_modificado.chr\n");
    fp = fopen("input_modificado.chr","w"); 
    fwrite(input,sizeof(unsigned char),max1*max2,fp);     
    fclose(fp);
  }

  if (debug ==1) {
    printf("Writing ouput laplacefile.flt\n");  
    fp = fopen("laplacefile.flt","w");
    for (i=0;i<max1;i++) {
      fwrite(laplacefield[i],sizeof(float),max2,fp);    
    }
    fclose(fp);
  }

  /* Reserve memory for maps and gradients*/
  maps = (float*)malloc(sizeof(float)*max1*max2);
  gradientx = (float**)malloc(sizeof(float*)*max1);
  for (i=0;i<max1;i++) {
    gradientx[i] = (float*)malloc(sizeof(float)*max2);
  }
  gradienty = (float**)malloc(sizeof(float*)*max1);
  for (i=0;i<max1;i++) {
    gradienty[i] = (float*)malloc(sizeof(float)*max2);
  }

  /*GradX(laplacefield, gradientx, max1, max2);
    GradY(laplacefield, gradienty, max1, max2);*/
  
  iGradX(laplacefield, gradientx, max1, max2);
  iGradY(laplacefield, gradienty, max1, max2);

  normalize(gradientx,gradienty,max1,max2);
  
  /* CODIGO DE CONTROL */
  if (debug == 1) {
    fp = fopen("gradientXnew.flt","w");
    if (fp == NULL) { 
      fprintf(stderr,"Failed writing gradientX.flt\n");
      exit(1); 
    } 
    for (i=0;i<max1;i++) {
      fwrite(gradientx[i],sizeof(float),max2,fp); 
    }
    fclose(fp); 
    
    fp = fopen("gradientYnew.flt","w"); 
    if (fp == NULL) { 
      fprintf(stderr,"Failed writing gradientY.flt\n");
      exit(1); 
    } 
    for (i=0;i<max1;i++) {
      fwrite(gradienty[i],sizeof(float),max2,fp);
    }
    fclose(fp);
  }
  /* CODIGO DE CONTROL */

  numasignaciones = 0;
  printf("Entering in thickness2D\n");
  gettimeofday(&endinit,NULL);
  if (reverse == 0) {
    if ( thickness2DYezzi(input, max1, max2, maps, laplacefield, gradientx, gradienty, num_it,hx,hy, label, debug) == 1 ) {
      printf("Error in thickness2D\n");
    }
  } else {
    if ( thickness2DYezzi_reverse(input, max1, max2, maps, laplacefield, gradientx, gradienty, num_it,hx,hy, label, debug) == 1 ) {
    printf("Error in thickness2D\n");
    }
  }
  /* if ( thickness2Dgradient(input, max1, max2, maps, laplacefield, gradientx, gradienty) == 1 ) {
    printf("Error in thickness2D\n");
    } */
  
  gettimeofday(&endtotal,NULL);

  printf("Writing ouput %s:\n",outputfile);
  write_float_output(outputfile, maps, max2, max1, 1, color_mode);

  free(laplacefield);
  free(maps);
  free(input);
  free(gradientx);
  free(gradienty);

  fprintf(stdout,"Initialization time: ");
  print_timing(stdout, startinit, endinit);
  fprintf(stdout,"thickness time: ");
  print_timing(stdout, endinit, endtotal);
  printf("OK thickness2D\n");
}
