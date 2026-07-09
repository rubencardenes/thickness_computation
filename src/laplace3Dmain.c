/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "laplace3D.h"

int main(int argc,char* argv[]) {
  unsigned char* input;
  float ***output;
  int i,j,k;
  int max1 = 256;
  int max2 = 256;
  int max3 = 70;
  int iterations = 10;
  float lambda = 0.5;
  FILE *fp,*fg;
  if (argc != 8 && argc != 7) {
    printf("Usage: laplace3D max1 max2 max3 input.vol output.volf iterations [lambda]\n");
    return 1;
  }

  max1 = atoi(argv[1]);
  max2 = atoi(argv[2]);
  max3 = atoi(argv[3]);
  iterations = atoi(argv[6]);
  if (argc == 8) {
    lambda = atof(argv[7]);
  }
  
  input = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2*max3);
  fp = fopen(argv[4],"rb");
  if (fp == NULL) {
    fprintf(stderr,"Failed reading inputfile %s\n",argv[4]);
    exit(1);
  }
  fread(input,sizeof(unsigned char),max1*max2*max3,fp);
  fclose(fp);
  
  output = (float***)malloc(sizeof(float**)*max3);
  for (k=0;k<max3;k++) {
    output[k] = (float**)malloc(sizeof(float)*max1);
  }
  for (k=0;k<max3;k++) {
    for (i=0;i<max1;i++) {
      output[k][i] = (float*)malloc(sizeof(float)*max2);
    }
  }
  
  printf("Entering in EdgeDetect3D\n");
  if ( EdgeDetect3D(input, max1, max2, max3) == 1 ) {
    printf("Error in EdgeDetect3D\n");
    return 1;
  } 
  
  printf("Relabeling\n");
  if ( relabel(input,max1*max2*max3,0,2) != 0) {
    printf("Error in Relabel\n");
    return 1;
  }
  
  printf("Entering in RelabelBoundary3D\n");
  if ( RelabelBoundary3D(input, max1, max2, max3) == 1 ) {
    printf("Error in RelabelBoundary\n");
    return 1;
  } 
  
  printf("Writing domain domain_modificado3d.vol\n");
  fg=fopen("../phantoms/thickness/domain_modificado3d.vol","w");
  fwrite(input,sizeof(unsigned char),max1*max2*max3,fg);    
  fclose(fg);

  printf("Entering in laplacian3D\n");
  if ( laplace3D(input, max1, max2, max3, output, iterations, lambda) == 1 ) {
    printf("Error in laplace3D\n");
    
  }

  printf("Writing ouput %s:\n",argv[5]);
  fg=fopen(argv[5],"w");
  for (k=0;k<max3;k++) {
    for (i=0;i<max1;i++) {
      fwrite(output[k][i],sizeof(float),max2,fg);    
    }
  }
  fclose(fg);

  free(input);
  free(output);
  printf("OK laplace2D\n");
  return 0;
}

