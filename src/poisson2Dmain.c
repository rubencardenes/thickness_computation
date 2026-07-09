/*(c) Ruben Cardenes Almeida, Valladolid, 8/11/2007 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <poisson2D.h>
#include <laplace2D.h>

int main(int argc,char* argv[]) {
  unsigned char* input,*min_local;
  float **output;
  int i;
  int max1 = 256;
  int max2 = 256;
  int iterations = 10;
  float lambda = 0.5,h = 1.0;
  FILE *fp,*fg;
  if (argc != 7 && argc != 6) {
    printf("Usage: poisson2D max1 max2 input.chr output.pnm iterations [lambda]\n");
    return 1;
  }

  int cortex_label = 2;

  max1 = atoi(argv[1]);
  max2 = atoi(argv[2]);
  iterations = atoi(argv[5]);
  if (argc == 7) {
    lambda = atof(argv[6]);
  }
  
  input = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2);
  fp = fopen(argv[3],"rb");
  if (fp == NULL) {
    fprintf(stderr,"Failed reading inputfile %s\n",argv[3]);
    exit(1);
  }
  fread(input,sizeof(unsigned char),max1*max2,fp);
  fclose(fp);
  
  output = (float**)malloc(sizeof(float*)*max1);
  for (i=0;i<max1;i++) {
    output[i] = (float*)malloc(sizeof(float)*max2);
  }
  
  printf("Entering in EdgeDetect\n");
  if ( EdgeDetect(input, max1, max2) == 1 ) {
    printf("Error in EdgeDetect\n");
  } 
  
  printf("Relabeling\n");
  if ( relabel(input,max1*max2,1,100) != 0) {
    printf("Error in Relabel\n");
  }
  
  printf("Writing domain domain_anillo_modificado.chr\n");
  fg=fopen("domain_anillo_modificado.chr","w");
  fwrite(input,sizeof(unsigned char),max1*max2,fg);    
  fclose(fg);

  printf("Entering in poisson2D\n");
  if ( poisson2D(input, max1, max2, output, iterations, lambda, 0, h, cortex_label) == 1 ) {
    printf("Error in poisson2D\n");
  }

  min_local = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2);
  printf("Entering in min local detection\n");
  if ( minimos_locales2D(output, min_local, max1, max2) == 1 ) {
    printf("Error in minimos_locales2D\n");
  }

  printf("Writing ouput %s:\n",argv[4]);
  fg=fopen(argv[4],"w");
  
  /* PGM */ 
  /* fprintf(fg,"P5\n");
  fprintf(fg,"# Generado por R. Cardenes\n");
  fprintf(fg,"%d %d\n",max1,max2);
  fprintf(fg,"255\n");*/
  
  /* MHD */
  /*fprintf(fg,"NDims = 2\n");
  fprintf(fg,"DimSize = %d %d\n",max1,max2);
  fprintf(fg,"ElementSize = 1 1\n");
  fprintf(fg,"ElementType = MET_FLOAT\n");
  fprintf(fg,"ElementByteOrderMSB = False\n");*/
  for (i=0;i<max1;i++) {
    fwrite(output[i],sizeof(float),max2,fg);    
  }
  fclose(fg);

  printf("Writing max_local.chr\n");
  fg=fopen("max_local.chr","w");
  fwrite(min_local,sizeof(unsigned char),max2*max2,fg);    
  fclose(fg);

  free(min_local);
  free(input);
  free(output);
  printf("OK poisson2D\n");
  return 0;
}

