/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include "io.h"
#include "thickness3D.h"
#include "laplace3D.h"
#include "DToptimo3d.h"
#include "png_write.h"

int numrechazos = 0;
int numasignaciones = 0;
int asignacionesraras = 0;
int numPrototypes;

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
	if (k > 0 && j > 0 && i > 0 && k < max3-1 && j < max2-1 && i < max1-1) {
	  if ( (segmented[sum]==label_wm) &&
	       ((segmented[sum+1]==label_cortex)||
		(segmented[sum-1]==label_cortex)||
		     
		(segmented[sum+max1]==label_cortex)||
		(segmented[sum-max1]==label_cortex)||
		     
		(segmented[sum+max1*max2]==label_cortex)||
		(segmented[sum-max1*max2]==label_cortex))) {
	      segmented[sum] = 1; 
	  } else if ( (segmented[sum]==0) &&
		      ((segmented[sum+1]==label_cortex)||
		       (segmented[sum-1]==label_cortex)||
		     
		       (segmented[sum+max1]==label_cortex)||
		       (segmented[sum-max1]==label_cortex)||
		     
		       (segmented[sum+max1*max2]==label_cortex)||
		       (segmented[sum-max1*max2]==label_cortex))) {	  
	    segmented[sum] = 128; 
	  }
	}
	sum++;
      }
    }
  }

  relabel_ushort(segmented,max1*max2*max3,0,255);
  relabel_ushort(segmented,max1*max2*max3,label_wm,255);
  relabel_ushort(segmented,max1*max2*max3,label_cortex,2);
  relabel_ushort(segmented,max1*max2*max3,128,0);
  return 0;
}

int main(int argc, char* argv[]) {
  unsigned char *input;
  unsigned short *input_short;
  float *maps = (float*)NULL,*maps_reverse = (float*)NULL;
  float **dist_maps = (float**)NULL;
  float **dist_maps_reverse = (float**)NULL;
  float ***laplacefield,***gradientx,***gradienty,***gradientz;
  float lambda = 0.5, mean, sigma =0;
  float hx =1,hy = 1,hz = 1;
  int i,j,k,col,row,c,option_index,swapbyte = 0, compute_mean = 0, thickness_DT = 0;
  int max1 = 256, max2 = 256, max3 = 1;
  int num_it = 10,iterations_laplace = 100, reverse = 0, suma = 0;
  int debug = 0,knee = 0, label_wm = 3, label_cortex = 2;
  FILE *fp,*fg;
  struct timeval startinit;
  struct timeval endinit;
  struct timeval endtotal;
  char input_prefix[200];
  char inputfile[200],outputfile[200],extension[6];
  char *laplacefile;
  unsigned char *aux,*l;

  while (1) {
    static struct option long_options[] = {
      {"hx", 1, 0, 0},
      {"hy", 1, 0, 0},
      {"hz", 1, 0, 0},
      {"lw", 1, 0, 0},
      {"lc", 1, 0, 0},
      {"DT", 0, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "dn:i:wrsl:km",long_options, &option_index);

    if (c == -1) {
      break;
    }               
      
    switch (c) {
    case 0:      
      if (strcmp(long_options[option_index].name,"hx") == 0) {
	hx = atof(optarg);
      }
      if (strcmp(long_options[option_index].name,"hy") == 0) {
	hy = atof(optarg);
      }
      if (strcmp(long_options[option_index].name,"hz") == 0) {
	hz = atof(optarg);
      }
      if (strcmp(long_options[option_index].name,"lc") == 0) {
	label_cortex = atoi(optarg);
      }
      if (strcmp(long_options[option_index].name,"lw") == 0) {
	label_wm = atoi(optarg);
      }
      if (strcmp(long_options[option_index].name,"DT") == 0) {
	thickness_DT = 1;
      }
      /* printf ("option %s = %f\n", long_options[option_index].name,threshold1);*/
      break;
    case 'd':
      debug = 1;   
      break;
    case 'n':
      num_it = atoi(optarg);
      break;
    case 'i':    
      iterations_laplace = atoi(optarg);
      break;    
    case 'w':
      swapbyte = 1;
      break;
    case 'r':
      reverse = 1;
      break;
    case 's':
      suma = 1;
      break;
    case 'l':
      lambda = atof(optarg);
      break;
    case 'k':
      knee = 1;
      break;
    case 'm':
      compute_mean = 1;
      break;
    case '?':
      printf("Author: Ruben Cardenes, April 2004\n");
      printf("Usage: thickness3D_cortex [options] segmented.vols/segmented_prefix output3D.volf nrows ncols nslices\n");
      printf("              -d (debug mode)\n");
      printf("              -n iterations_thickness (10)\n");
      printf("              -i iterations_laplace (100)\n");
      printf("              -w (swapbytes) \n");
      printf("              -r (reverse) \n");
      printf("              -s (sum) \n");
      printf("              -k (knee mode) \n");
      printf("              -l lambda (0.5)\n");
      printf("              -m (compute and show mean thickness)\n");
      printf("              --lw white matter label \n");
      printf("              --lc cortex label");
      printf("              --hx hy (1) \n");
      printf("              --hy hx (1) \n");
      printf("              --hz hz (1) \n");
      printf("              --DT (compute thickness with Euclidean DT)\n");
      return 1;
      break;
    
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  
  if ((argc - optind) != 5) {
    printf ("Incorrect number of arguments: ");
    printf("Author: Ruben Cardenes, April 2004 \n");
    printf("Usage: thickness3D_cortex  [options] segmented.vols/segmented_prefix output3D.volf nrows ncols nslices\n");
    printf("              -d (debug mode)\n");
    printf("              -n iterations_thickness (10)\n");
    printf("              -i iterations_laplace (100)\n");
    printf("              -w (swapbytes) \n");
    printf("              -r (reverse) \n");
    printf("              -s (sum) \n");
    printf("              -k (knee mode) \n");
    printf("              -l lambda (0.5)\n");
    printf("              -m compute and show mean thickness\n");
    printf("              --lw white matter label \n");
    printf("              --lc cortex label\n");
    printf("              --hx hy (1) \n");
    printf("              --hy hx (1) \n");
    printf("              --hz hz (1) \n");
    printf("              --DT (compute thickness with Euclidean DT)\n");
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
  input = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2*max3);
  input_short = (unsigned short*)malloc(sizeof(unsigned short)*max1*max2*max3);

  for (i=0;i<max1*max2*max3;i++) {
    input[i] = -1;
  }

  l = strchr(input_prefix,'.');
  if (l != (unsigned char*)NULL) {
    strncpy(extension,(char*)l,5);
    extension[5] = '\0';
  } else {
    strcpy(extension,"");
  }
  if (strncmp(extension,".vols",5) == 0) {
    fp = fopen(input_prefix,"r");
    if (fp == NULL) {
      fprintf(stderr,"Failed reading inputfile %s\n",input_prefix);
      exit(1);
    }
    fread(input_short,sizeof(unsigned short),max1*max2*max3,fp);
    fclose(fp);
  } else {
    /* Read data */
    for (i=0;i<max3;i++) {
      sprintf(inputfile,"%s.%03d",input_prefix,i+1);
      fp = fopen(inputfile,"r");
      if (fp == NULL) {
	fprintf(stderr,"Failed reading inputfile %s\n",inputfile);
	exit(1);
      }
      fread(&input_short[max1*max2*i],sizeof(unsigned short),max1*max2,fp);
      fclose(fp);
    }
  }
  
  if (swapbyte == 1) {
    for (i=0; i < max1*max2*max3; i++) {
      aux = (unsigned char*)&input_short[i];
      input_short[i] = ReadGEShort(aux);
    }
  }

  /* Report the distinct label values in the volume. In cortex mode require the
     given --lw and --lc to be present, so a wrong label fails loudly instead of
     silently producing an empty result. In knee mode the boundaries are already
     encoded in the input, so --lw/--lc are not validated. */
  {
    unsigned char present[256];
    print_domain_values_ushort(input_short, max1*max2*max3, present);
    if (knee == 0) {
      if (!require_label(present, label_wm, "--lw") ||
          !require_label(present, label_cortex, "--lc")) {
        exit(1);
      }
    }
  }

  if (knee == 1) {
    EdgeDetect3D_knee(input, max1, max2, max3);
  } else {
    compute_boundary_cortex3D(input_short,max1, max2, max3, label_wm, label_cortex);
  }

  /* After boundary processing the band is labeled 2 (cortex mode relabels the
     input --lc to 2; knee mode expects the band to already be 2). The thickness
     computation uses that fixed label internally, so pin label_cortex to 2 here
     as well: this keeps the -m mean statistics correct regardless of the input
     --lw/--lc values. */
  label_cortex = 2;

  for (i=0;i<max1*max2*max3;i++) {
    input[i] = (unsigned char)input_short[i];
  }
  free(input_short);

  if (thickness_DT == 0) {

    laplacefield = (float***)malloc(sizeof(float**)*max3);
    for (k=0;k<max3;k++) {
      laplacefield[k] = (float**)malloc(sizeof(float*)*max1);
    }
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	laplacefield[k][i] = (float*)malloc(sizeof(float)*max2);
      }
    }

    printf("Entering in laplacian3D\n");
    if ( laplace3D_voxelsize(input, max1, max2, max3, laplacefield, iterations_laplace, hx, hy, hz, lambda) == 1 ) {
      printf("Error in thickness3D\n");
    }
    /*if ( laplace3D(input, max1, max2, max3, laplacefield, iterations_laplace, lambda) == 1 ) {
      printf("Error in thickness3D\n");
      }*/

    if (debug == 1) {
      printf("Writing domain domain_modificado3d.vol\n");
      fp=fopen("domain_modificado3d.vol","w");
      fwrite(input,sizeof(unsigned char),max1*max2*max3,fp);    
      fclose(fp);

      printf("Writing ouput laplacefile.volf\n");  
      fp = fopen("laplacefile.volf","w");
      for (k=0;k<max3;k++) {
	for (i=0;i<max1;i++) {
	  fwrite(laplacefield[k][i],sizeof(float),max2,fp);    
	}
      }
      fclose(fp);       
    }

    /* Reserve memory for maps and gradients*/
    maps = (float*)malloc(sizeof(float)*max1*max2*max3);
    gradientx = (float***)malloc(sizeof(float**)*max3);
    for (k=0;k<max3;k++) {
      gradientx[k] = (float**)malloc(sizeof(float*)*max1);
    }
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	gradientx[k][i] = (float*)malloc(sizeof(float)*max2);
      }
    }
    
    gradienty = (float***)malloc(sizeof(float**)*max3);
    for (k=0;k<max3;k++) {
      gradienty[k] = (float**)malloc(sizeof(float*)*max1);
    }
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	gradienty[k][i] = (float*)malloc(sizeof(float)*max2);
      }
    }

    gradientz = (float***)malloc(sizeof(float**)*max3);
    for (k=0;k<max3;k++) {
      gradientz[k] = (float**)malloc(sizeof(float*)*max1);
    }
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	gradientz[k][i] = (float*)malloc(sizeof(float)*max2);
      }
    }
    
    printf("Doing gradients \n");
    iGradX3D(laplacefield, gradientx, max2, max1, max3, hx);
    iGradY3D(laplacefield, gradienty, max2, max1, max3, hy);
    iGradZ3D(laplacefield, gradientz, max2, max1, max3, hz);

    normalize3D(gradientx,gradienty,gradientz,max1,max2,max3);
  
    /* CODIGO DE CONTROL */
    if (debug == 1) {
      fp = fopen("gradientX.volf","w");
      if (fp == NULL) { 
	fprintf(stderr,"Failed writing gradientX.volf\n");
	exit(1); 
      } 
      for (k=0;k<max3;k++) {
	for (i=0;i<max1;i++) {
	  fwrite(gradientx[k][i],sizeof(float),max2,fp);    
	}
      }
      fclose(fp); 
      
      fp = fopen("gradientY.volf","w"); 
      if (fp == NULL) { 
	fprintf(stderr,"Failed writing gradientY.volf\n");
	exit(1); 
      } 
      for (k=0;k<max3;k++) {
	for (i=0;i<max1;i++) {
	  fwrite(gradienty[k][i],sizeof(float),max2,fp);    
	}
      }
      fclose(fp);

      fp = fopen("gradientZ.volf","w"); 
      if (fp == NULL) { 
	fprintf(stderr,"Failed writing gradientZ.volf\n");
	exit(1); 
      } 
      for (k=0;k<max3;k++) {
	for (i=0;i<max1;i++) {
	  fwrite(gradientz[k][i],sizeof(float),max2,fp);    
	}
      }
      fclose(fp);
    }
    /* CODIGO DE CONTROL */

    printf("Entering in thickness3D\n");
    gettimeofday(&endinit,NULL);
    if (suma == 0) {
      if (reverse == 0) {
	if ( thickness3DYezzi(input, max1, max2, max3, maps, laplacefield, gradientx, gradienty, gradientz, num_it,hx,hy,hz) == 1 ) {
	  printf("Error in thickness3D\n");
	}
      } else {
	if ( thickness3DYezzi_reverse(input, max1, max2, max3,  maps, laplacefield, gradientx, gradienty, gradientz, num_it,hx,hy,hz) == 1 ) {
	  printf("Error in thickness3D\n");
	}
      }
    }

    if (suma == 1) {
      maps_reverse = (float*)malloc(sizeof(float)*max1*max2*max3);
      if ( thickness3DYezzi(input, max1, max2, max3, maps, laplacefield, gradientx, gradienty, gradientz, num_it,hx,hy,hz) == 1 ) {
	printf("Error in thickness3D\n");
      }
      if ( thickness3DYezzi_reverse(input, max1, max2, max3,  maps_reverse, laplacefield, gradientx, gradienty, gradientz, num_it,hx,hy,hz) == 1 ) {
	printf("Error in thickness3D\n");
      }
      relabel_float(maps,max1*max2*max3,-1,0);
      relabel_float(maps_reverse,max1*max2*max3,-1,0);
      sumar_l1l2(maps,maps_reverse,maps,max1*max2*max3);     
      free(maps_reverse);
    }

    relabel_float(maps,max1*max2*max3,-1,0);
    printf("Writing output %s:\n",outputfile);
    fp=fopen(outputfile,"w"); 
    if (fp == NULL) { 
      fprintf(stderr,"Failed writing %s\n",outputfile);
      return 1; 
    } 
    fwrite(maps,sizeof(float),max1*max2*max3,fp);    
    fclose(fp);
    
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	free(laplacefield[k][i]);
	free(gradientx[k][i]);
	free(gradienty[k][i]);
	free(gradientz[k][i]);
      }
      free(laplacefield[k]);
      free(gradientx[k]);
      free(gradienty[k]);
      free(gradientz[k]);
    }
    free(laplacefield);
    free(gradientx);
    free(gradienty);
    free(gradientz);
  } else {
    dist_maps = (float**)malloc(sizeof(float*)*3);
    dist_maps[0] = (float*)malloc(sizeof(float)*max1*max2*max3);    
    dist_maps[1] = (float*)malloc(sizeof(float)*max1*max2*max3);
    dist_maps[2] = (float*)malloc(sizeof(float)*max1*max2*max3);
   
    if (suma == 0) {
      if (reverse == 1) {
	DToptimo3d(input, max1, max2, max3, 1, dist_maps, 1, 2, debug, hx,hy,hz);
      } else {
	DToptimo3d(input, max1, max2, max3, 1, dist_maps, 0, 2, debug, hx,hy,hz);
      }
      relabel_float(dist_maps[2],max1*max2*max3,999999,0);
    }
    if ( suma == 1) { 
       dist_maps_reverse = (float**)malloc(sizeof(float*)*3);
       dist_maps_reverse[0] = (float*)malloc(sizeof(float)*max1*max2*max3);    
       dist_maps_reverse[1] = (float*)malloc(sizeof(float)*max1*max2*max3);
       dist_maps_reverse[2] = (float*)malloc(sizeof(float)*max1*max2*max3);
       DToptimo3d(input, max1, max2, max3, 1, dist_maps, 1, 2, debug, hx,hy,hz);
       DToptimo3d(input, max1, max2, max3, 1, dist_maps_reverse, 0, 2, debug, hx,hy,hz);
       relabel_float(dist_maps[2],max1*max2*max3,999999,0);
       relabel_float(dist_maps_reverse[2],max1*max2*max3,999999,0);
       sumar_l1l2(dist_maps[2],dist_maps_reverse[2],dist_maps[2],max1*max2*max3); 
       free(dist_maps_reverse);
    }

    printf("Writing output %s:\n",outputfile);
    fp=fopen(outputfile,"w"); 
    if (fp == NULL) { 
      fprintf(stderr,"Failed writing %s\n",outputfile);
      return 1; 
    }
    fwrite(dist_maps[2],sizeof(float),max1*max2*max3,fp);    
    fclose(fp);
  }
  
  gettimeofday(&endtotal,NULL);
  if (compute_mean) {
    if (suma == 0) {
      if (thickness_DT == 0) {
	if (reverse) { 
	  mean = compute_mean_thickness(input,maps,label_cortex,0,max1,max2,max3,&sigma);
	} else {
	  mean = compute_mean_thickness(input,maps,label_cortex,1,max1,max2,max3,&sigma);
	} 
      } else {
	if (reverse) { 
	  mean = compute_mean_thickness(input,dist_maps[2],label_cortex,0,max1,max2,max3,&sigma);
	} else {
	  mean = compute_mean_thickness(input,dist_maps[2],label_cortex,1,max1,max2,max3,&sigma);
	}
      }
      printf("Mean thickness in the surface boundary = %f std = %f\n",mean,sigma);
    } else {
      if (thickness_DT == 0) {
	mean = compute_mean_thickness_volume(input,maps,label_cortex,max1,max2,max3,&sigma); 
      }	else {
	mean = compute_mean_thickness_volume(input,dist_maps[2],label_cortex,max1,max2,max3,&sigma); 
      }
      printf("Mean thickness in the volume boundary = %f std = %f\n",mean,sigma);
      
    }
  }

  if (maps != (float*)NULL) {
    free(maps);
  }
  if (dist_maps != (float**)NULL) {
    free(dist_maps);
  }
  free(input);
  
  fprintf(stdout,"Initialization ");
  print_timing(stdout, startinit, endinit);
  fprintf(stdout,"thickness ");
  print_timing(stdout, endinit, endtotal);
  printf("OK thickness3D\n");
  return 0;
}
