/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include "io.h"
#include "laplace2D.h"

int compute_boundary_cortex3D(unsigned short *segmented,int max1,int max2, int max3, int label_wm, int label_cortex) {
  int i,j,k,sum,found;

  for (i=0;i<max1*max2*max3;i++) {
    if (segmented[i] != label_cortex && segmented[i] != label_wm && segmented[i] != 0 ) {
      relabel_ushort(segmented,max1*max2*max3,segmented[i],0);
    }
  }

  sum=0;
  for(k=0; k<max3; k++) {
    for(j=0; j<max2; j++) {
      for(i=0; i<max1; i++) {     
	if ((i==0)||(j==0)||(k==0)||(i==max1-1)||(j==max2-1)||(k==max3-1) ) {
	  /* nothing to do */	 
	} else if ( (segmented[sum]==label_wm) &&
		    ((segmented[sum+1]==label_cortex)||
		     (segmented[sum-1]==label_cortex)||
		     
		     (segmented[sum+max2]==label_cortex)||
		     (segmented[sum-max2]==label_cortex)||
		     
		     (segmented[sum+max1*max2]==label_cortex)||
		     (segmented[sum-max1*max2]==label_cortex))) {
	  segmented[sum] = 1; 
	} else if ( (segmented[sum]==0) &&
		    ((segmented[sum+1]==label_cortex)||
		     (segmented[sum-1]==label_cortex)||
		     
		     (segmented[sum+max2]==label_cortex)||
		     (segmented[sum-max2]==label_cortex)||
		     
		     (segmented[sum+max1*max2]==label_cortex)||
		     (segmented[sum-max1*max2]==label_cortex))) {	  
	  segmented[sum] = 128; 
	}
	sum++;
      }
    }
  }

  return 0;
}

int main(int argc, char* argv[]) {
  unsigned short *input;
  int i,j,swapbyte,label_wm,label_cortex,tam,fsize,hsize,c,option_index;
  float threshold;
  int max1 = 256;
  int max2 = 256;
  int max3 = 1;
  int debug = 1;
  FILE *fp,*fg;
  struct timeval startinit;
  struct timeval endinit;
  struct timeval endtotal;
  char input_prefix[200];
  char input_file[200],outputfile[200],extension[5];
  unsigned char *aux,*output,*l;
  float *output_float;

  while (1) {
    static struct option long_options[] = {
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "dwl:c:",long_options, &option_index);

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
      label_wm = atoi(optarg);
      break;
    case 'c':
      label_cortex = atoi(optarg);
      break; 
    case '?':
      printf("Author: Ruben Cardenes \n");
      printf("Usage: compute_boundary_cortex3D [options] input_prefix output3D.vols nrows ncols nslices\n");
      printf("              -d (debug mode)\n");
      printf("              -w (swapbytes) \n");
      printf("              -l label_wm \n");
      printf("              -c label_cortex \n");
   
      return 1;
      break;
    
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  
  if ((argc - optind) != 5) {
    printf ("Incorrect number of arguments: ");
    printf("Author: Ruben Cardenes \n");
    printf("Usage: compute_boundary_cortex3D [options] input_prefix output3D.vols nrows ncols nslices\n");
    printf("              -d (debug mode)\n");
    printf("              -w (swapbytes) \n");
    printf("              -l label_wm \n");
    printf("              -c label_cortex \n");   
    return 1;
  } else {
    while (optind < argc) {
      if (sscanf(argv[optind++], "%s", input_prefix) == 0)
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
  input = (unsigned short*)malloc(sizeof(unsigned short)*max1*max2*max3);

  l = strchr(input_prefix,'.');
  strncpy(extension,l,5);
  if (strcmp(extension,".vols") == 0) {
    fp = fopen(input_prefix,"r");
    if (fp == NULL) {
      fprintf(stderr,"Failed reading inputfile %s\n",input_prefix);
      exit(1);
    }
    fread(input,sizeof(unsigned short),max1*max2*max3,fp);
    fclose(fp);
  } else {
    /* Read data */
    for (i=0;i<max3;i++) {
      sprintf(input_file,"%s.%03d",input_prefix,i+1);
      fp = fopen(input_file,"r");
      if (fp == NULL) {
	fprintf(stderr,"Failed reading inputfile %s\n",input_file);
	exit(1);
      }
      fread(&input[max1*max2*i],sizeof(unsigned short),max1*max2,fp);
      fclose(fp);
    }
  }

  if (swapbyte == 1) {
    for (i=0; i < max1*max2*max3; i++) {      
      aux = (unsigned char*)&input[i];
      input[i] = ReadGEShort(aux);   
    }
  }
  
  /* Compute */
  compute_boundary_cortex3D(input,max1,max2,max3,label_wm,label_cortex);

  relabel_ushort(input,max1*max2*max3,0,255);
  relabel_ushort(input,max1*max2*max3,label_wm,255);
  relabel_ushort(input,max1*max2*max3,label_cortex,2);
  relabel_ushort(input,max1*max2*max3,128,0);
  /* Write data */
  printf("writing file %s\n",outputfile);
  fp = fopen(outputfile,"w");
  if (fp == NULL) {
    fprintf(stderr,"Failed writing file %s\n",outputfile);
    return 1;
  }
  fwrite(input,sizeof(unsigned short),max1*max2*max3,fp);  
  fclose(fp);

  /* Free data */
  free(input);
  printf("compute_boundary_cortex3D OK\n");
  return 0;

}
