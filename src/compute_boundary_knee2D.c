/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "io.h"
#include "laplace2D.h"

int test_boundary(unsigned short *segmented,unsigned short *original,int i,int j,int max2,float threshold) {
  int x,y,mapindex;
  int count = 0;
  int sum = 0;
  float media;
  int n_size = 6;
  mapindex = i*max2+j;
  for(x=-n_size; x<=n_size; x++) {
    for(y=-n_size; y<=n_size; y++) {
      if (x ==0 && y==0) continue;
      if (segmented[mapindex+max2*y+x] != segmented[mapindex]) {
	count++;
	sum += original[mapindex+max2*y+x];
      }
    }
  }

  media = (float)sum/(float)count;
  if (media < threshold) {
    segmented[mapindex] = 1; 
  }
  else { 
    segmented[mapindex] = 2; 
  }
 
  return 0;
}

int compute_boundary2D(unsigned short *segmented,unsigned short* original,int max1,int max2,int label, float threshold) {
  int i,j,istart,jstart,sum;
  
  istart = -1;
  jstart = -1;
  sum = 0;
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (segmented[sum] == label) { 
	istart = i; jstart =j;	
      }
      if (segmented[sum] != 0 && segmented[sum] != label) { 	
	relabel_ushort(segmented,max1*max2,segmented[sum],0);
      }
      sum++;
    }
  }
  
  if (istart == -1 || jstart == -1) {
    printf("Label %d NOT found \n",label); 
    return 0;
  }

  printf("label = %d found istart %d jstart %d \n",label,istart,jstart);
  sum=0;
  for(i=0; i<max1; i++) {
    for(j=0; j<max2; j++) {
      if ((i==0)||(j==0)||(i==max1-1)||(j==max2-1)) {
	/* nothing to do */	 
      } else if ( (segmented[sum]==label)&&
		  ((segmented[sum+1]==0)||
		   (segmented[sum-1]==0)||
		   
		   (segmented[sum+max2]==0)||
		   (segmented[sum-max2]==0))) {	   
	   test_boundary(segmented,original,i,j,max2,threshold);
      }
      sum++;
    }
  }
  
  return 0;
}

int main(int argc, char* argv[]) {
  unsigned short *segmented;
  unsigned short *original;
  int i,j,swapbyte,label,tam,fsize,hsize;
  float threshold;
  int max1 = 256;
  int max2 = 256;
  int max3 = 1;
  int debug = 1;
  FILE *fp,*fg;
  struct timeval startinit;
  struct timeval endinit;
  struct timeval endtotal;
  char *segmented_file,*original_file;
  unsigned char *aux,*output;
  float *output_float;

  if (argc != 9) {
    printf("Usage: compute_boundary2D max1 max2 segmented_file.ush original_file.mri output2D.chr label threshold swapbyte(0/1) \n");
    return 1;
  }

  max1 = atoi(argv[1]);
  max2 = atoi(argv[2]);
  segmented_file = argv[3];
  original_file = argv[4];
  label = atoi(argv[6]);
  threshold = atof(argv[7]);
  swapbyte = atoi(argv[8]);

  gettimeofday(&startinit,NULL);
  /* reserve data */
  segmented = (unsigned short*)malloc(sizeof(unsigned short)*max1*max2);  
  original = (unsigned short*)malloc(sizeof(unsigned short)*max1*max2);

  /* Read data */
  fp = fopen(segmented_file,"r");
  if (fp == NULL) {
    fprintf(stderr,"Failed reading inputfile %s\n",segmented_file);
    exit(1);
  }
  fread(segmented,sizeof(unsigned short),max1*max2,fp);
  fclose(fp);

  fsize = fileSize(original_file);
  hsize = headerSize(fsize);

  fp = fopen(original_file,"r");
  if (fp == NULL) {
    fprintf(stderr,"Failed reading inputfile %s\n",original_file);
    exit(1);
  }
  if (hsize >0) {
    fseek(fp, hsize, 1);
  }
  fread(original,sizeof(unsigned short),max1*max2,fp);
  fclose(fp);
  /* If we have to swap the data: */
  if (swapbyte == 1) {
    for (i=0; i < max1*max2; i++) {
	aux = (unsigned char*)&segmented[i];
	segmented[i] = ReadGEShort(aux);      
	aux = (unsigned char*)&original[i];
	original[i] = ReadGEShort(aux);   
    }
  }

  /* Compute */
  compute_boundary2D(segmented,original,max1,max2,label,threshold);
  sizefilter2D(segmented, max1,max2, 15, 1, 2);
  sizefilter2D(segmented, max1,max2, 20, 2, 1); 
  tam = maxcomponent2D(segmented, max1, max2, 1); 
  sizefilter2D(segmented, max1,max2, tam-1, 1, 2); 

  relabel_ushort(segmented, max1*max2, 0, 255);
  relabel_ushort(segmented, max1*max2, 2, 0);
  relabel_ushort(segmented, max1*max2, label, 2);
  
  printf("computing corners\n");
  new_compute_corners(segmented, max1,max2);

  /* convert data */
  output=(unsigned char*)malloc(max1*max2*sizeof(unsigned char));
  for (i=0;i<max1*max2;i++) {
    output[i] = segmented[i];
  }

  /* control code */
  output_float=(float*)malloc(max1*max2*sizeof(float));
  for (i=0;i<max1*max2;i++) {
    output_float[i] = segmented[i];
  }
  fp = fopen("/home/ruben/data/knee-images/thickness/010/output2D.flt","w");
  if (fp == NULL) {
    fprintf(stderr,"Failed writing file \n");
    exit(1);
  }
  fwrite(output_float,sizeof(float),max1*max2,fp);
  fclose(fp);
  /* control code */

  /* Write data */
  fp = fopen(argv[5],"w");
  if (fp == NULL) {
    fprintf(stderr,"Failed writing file %s\n",argv[5]);
    exit(1);
  }
  fwrite(output,sizeof(unsigned char),max1*max2,fp);
  
  fclose(fp);

  /* Free data */
  free(original); 
  free(segmented);
  free(output);
  free(output_float);
  printf("compute_boundary2D OK\n");
  return 0;

}
