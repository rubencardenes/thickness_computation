/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include "io.h"
#include "laplace2D.h"

int test_boundary3D(unsigned char *segmented,unsigned short *original,int mapindex,int max1,int max2,float threshold) {
  int x,y,z;
  int count0 = 0,count1 = 0;
  int sum = 0;
  float media;
  int n_size = 2;

  for(x=-1; x<=1; x++) {
    for(y=-1; y<=1; y++) {
      for (z=-1; z<=1; z++) {
	if (x==0 && y==0 && z==0) continue;
	if (segmented[mapindex+max1*max2*z+max1*y+x] == 1) {
	  count1++;
	}
	if (segmented[mapindex+max1*max2*z+max1*y+x] == 0) {
	  count0++;
	}
      }
    }
  }

  if ((count1+count0) != 0) {
    if (count1 > count0) {
      segmented[mapindex] = 1;
    } else {
      segmented[mapindex] = 0;
    }
  } else {
    
    for(x=-n_size; x<=n_size; x++) {
      for(y=-n_size; y<=n_size; y++) {
	for (z=-n_size; z<=n_size; z++) {
	  if (x==0 && y==0 && z==0) continue;
	  if (segmented[mapindex+max1*max2*z+max1*y+x] != segmented[mapindex]) {
	    count0++;
	    sum += original[mapindex+max1*max2*z+max1*y+x];
	  }
	}
      }
    }
    
    media = (float)sum/(float)count0;
    if (media < threshold) {
      segmented[mapindex] = 1; 
    }
    else { 
      segmented[mapindex] = 2; 
    }
  }
 
  return 0;
}

int compute_boundary3D(unsigned char *segmented,unsigned short* original,int max1,int max2, int max3, int label, float threshold) {
  int i,j,k,sum,found;
  
  /*found = 0;
  for (i=0;i<max1*max2*max3;i++) {
    if (segmented[i] == label) { 	
      found = 1;
    }
    if (segmented[i] != 0 && segmented[i] != label) { 	
      relabel(segmented,max1*max2*max3,segmented[i],0);
    }
  }

  if (found ==1) {
    printf("Found label %d \n",label);
    }*/

  sum=0;
  for(k=0; k<max3; k++) {
    for(j=0; j<max2; j++) {
      for(i=0; i<max1; i++) {     
	/* test_boundary3D's fallback path scans a window of radius n_size=2,
	   so the margin here must be at least 2 to keep that scan in bounds. */
	if ((i<2)||(j<2)||(k<2)||(i>=max1-2)||(j>=max2-2)||(k>=max3-2) ) {
	/* nothing to do */
	} else if ( (segmented[sum]==0 || segmented[sum]==1
		     || segmented[sum]==2)&&
		    ((segmented[sum+1]==255)||
		     (segmented[sum-1]==255)||
		     
		     (segmented[sum+max1]==255)||
		     (segmented[sum-max1]==255)||
		     
		     (segmented[sum+max1*max2]==255)||
		     (segmented[sum-max1*max2]==255))) {
	  test_boundary3D(segmented,original,sum,max1,max2,threshold);
	}
	sum++;
      }
    }
  }
  
  return 0;
}

int main(int argc, char* argv[]) {
  unsigned char *input;
  unsigned short *original;
  int i,j,swapbyte=0,label=0,tam,fsize,hsize,c,option_index;
  float threshold=0;
  int max1 = 256;
  int max2 = 256;
  int max3 = 1;
  int debug = 1;
  FILE *fp,*fg;
  struct timeval startinit;
  struct timeval endinit;
  struct timeval endtotal;
  char segmented_prefix[200],original_prefix[200];
  char segmented_file[200],original_file[200],outputfile[200];
  unsigned char *aux,*output;
  float *output_float;

  while (1) {
    static struct option long_options[] = {
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "dwl:t:",long_options, &option_index);

    if (c == -1) {
      break;
    }               
      
    switch (c) {
    case 0:      
     
      /* printf ("option %s = %f\n", long_options[option_index].name,threshold1);*/
      break;
    case 'd':
      debug = 1;   
      break;
    case 'w':
      swapbyte = 1;
      break;      
    case 'l':
      label = atoi(optarg);
      break;
    case 't':
      threshold = atof(optarg);
      break; 
    case '?':
      printf("Author: Ruben Cardenes \n");
      printf("Usage: compute_boundary_knee3D [options] segmented_prefix original_prefix output3D.vols nrows ncols nslices\n");
      printf("              -d (debug mode)\n");
      printf("              -w (swapbytes) \n");
      printf("              -l label \n");
      printf("              -t threshold \n");
      return 1;
      break;
    
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  
  if ((argc - optind) != 6) {
    printf ("Incorrect number of arguments: ");
    printf("Author: Ruben Cardenes \n");
    printf("Usage: compute_boundary_knee3D [options] segmented_prefix original_prefix output3D.vols nrows ncols nslices\n");
    printf("              -d (debug mode)\n");
    printf("              -w (swapbytes) \n");
    printf("              -l label\n");
    printf("              -t threshold \n");
    return 1;
  } else {
    while (optind < argc) {
      if (sscanf(argv[optind++], "%s", segmented_prefix) == 0)
	printf ("Error parsing argument \n");
      if (sscanf(argv[optind++], "%s", original_prefix) == 0)
	printf ("Error parsing argument \n");
      if (sscanf(argv[optind++], "%s", outputfile) == 0)
	printf ("Error parsing argument \n");
      if (sscanf(argv[optind++], "%d", &max1) == 0)
	printf ("Error parsing argument \n");
      if (sscanf(argv[optind++], "%d", &max2) == 0)
	printf ("Error parsing argument \n");
      if (sscanf(argv[optind++], "%d", &max3) == 0)
	printf ("Error parsing argument \n");      
    }
  }

  gettimeofday(&startinit,NULL);
  /* reserve data */
  input = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2*max3);
  original = (unsigned short*)malloc(sizeof(unsigned short)*max1*max2*max3);

  /* Read data */
  for (i=0;i<max3;i++) {
    sprintf(segmented_file,"%s.%03d",segmented_prefix,i+1);
    fp = fopen(segmented_file,"r");
    if (fp == NULL) {
      fprintf(stderr,"Failed reading inputfile %s\n",segmented_file);
      exit(1);
    }
    fread(&input[max1*max2*i],sizeof(unsigned char),max1*max2,fp);
    fclose(fp);
  }

  for (i=0;i<max3;i++) {
    sprintf(original_file,"%s.%03d",original_prefix,i+1);
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
    fread(&original[max1*max2*i],sizeof(unsigned short),max1*max2,fp);
    fclose(fp);
  }
  /* If we have to swap the data: */
  if (swapbyte == 1) {
    for (i=0; i < max1*max2*max3; i++) {      
      aux = (unsigned char*)&original[i];
      original[i] = ReadGEShort(aux);   
    }
  }

  /* Compute */
  compute_boundary3D(input,original,max1,max2,max3,label,threshold);
  /* sizefilter2D(segmented, max1,max2, 15, 1, 2);
  sizefilter2D(segmented, max1,max2, 20, 2, 1); 
  tam = maxcomponent2D(segmented, max1, max2, 1); 
  sizefilter2D(segmented, max1,max2, tam-1, 1, 2); 

  relabel_ushort(segmented, max1*max2, 0, 255);
  relabel_ushort(segmented, max1*max2, 2, 0);
  relabel_ushort(segmented, max1*max2, label, 2);*/
   
  /* control code */
  /* output_float=(float*)malloc(max1*max2*sizeof(float));
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
  free(output_float); */
  /* control code */

  /* Write data */
  printf("writing file %s\n",outputfile);
  fp = fopen(outputfile,"w");
  if (fp == NULL) {
    fprintf(stderr,"Failed writing file %s\n",outputfile);
    exit(1);
  }
  fwrite(input,sizeof(unsigned char),max1*max2*max3,fp);  
  fclose(fp);

  /* Free data */
  free(original); 
  free(input);
  printf("compute_boundary3D OK\n");
  return 0;

}
