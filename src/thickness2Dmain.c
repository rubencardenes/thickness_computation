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

int compute_boundary_cortex2D(unsigned char *segmented,int max1,int max2, int label_wm, int label_cortex) {
  int i,j,k,sum,found;

  for (i=0;i<max1*max2;i++) {
    if (segmented[i] != label_cortex && segmented[i] != label_wm && segmented[i] != 0 ) {
      relabel(segmented,max1*max2,segmented[i],0);
    }
  }

  sum=0;
  for(i=0; i<max1; i++) {
    for(j=0; j<max2; j++) {
      if ((i==0)||(j==0)||(i==max1-1)||(j==max2-1)) {
	  /* nothing to do */
	} else if ( (segmented[sum]==label_wm) &&
		    ((segmented[sum+1]==label_cortex)||
		     (segmented[sum-1]==label_cortex)||

		     (segmented[sum+max2]==label_cortex)||
		     (segmented[sum-max2]==label_cortex))) {
	  segmented[sum] = 1;
	} else if ( (segmented[sum]==0) &&
		    ((segmented[sum+1]==label_cortex)||
		     (segmented[sum-1]==label_cortex)||

		     (segmented[sum+max2]==label_cortex)||
		     (segmented[sum-max2]==label_cortex))) {
	  segmented[sum] = 128;
	}
      sum++;
    }
  }

  relabel(segmented,max1*max2,0,255);
  relabel(segmented,max1*max2,label_wm,255);
  relabel(segmented,max1*max2,label_cortex,2);
  relabel(segmented,max1*max2,128,0);
  return 0;
}

int main(int argc, char* argv[]) {
  unsigned char *input;
  float *maps,*maps_reverse,*streaml_maps;
  float **laplacefield,**gradientx,**gradienty;
  float lambda = 0.5, hx = 1,hy = 1;
  int i,col,row,c,option_index,num_it = 10,iterations_laplace = 100, suma = 0, swapbyte = 0;
  int max1 = 256;
  int max2 = 256;
  int max3 = 1;
  int reverse = 0,compute_mean = 0, label_cortex = 2, label_wm = 3;
  int debug = 0, streamlines = 0;
  int thickness_DT = 0;
  int color_mode = COLOR_GRAY;
  int width, height;
  FILE *fp,*fg;
  struct timeval startinit;
  struct timeval endinit;
  struct timeval endtotal;
  char inputfile[200],outputfile[200],*laplacefile;
  
  while (1) {
    static struct option long_options[] = {
      {"hx", 1, 0, 0},
      {"hy", 1, 0, 0},
      {"lw", 1, 0, 0},
      {"lc", 1, 0, 0},
      {"DT", 0, 0, 0},
      {"streamlines", 0, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "dn:i:wrsl:mc:",long_options, &option_index);

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
      if (strcmp(long_options[option_index].name,"lc") == 0) {
	label_cortex = atoi(optarg);
      }
      if (strcmp(long_options[option_index].name,"lw") == 0) {
	label_wm = atoi(optarg);
      }
      if (strcmp(long_options[option_index].name,"DT") == 0) {
	thickness_DT = 1;
      }
      if (strcmp(long_options[option_index].name,"streamlines") == 0) {
	streamlines = 1;
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
    case 'm':
      compute_mean = 1;
      break;
    case 'c':
      color_mode = atoi(optarg);
      break;
    case '?':
      printf("Author: Ruben Cardenes, April 2004\n");
      printf("Usage: thickness2D [options] input2D.png output2D.png\n");
      printf("              -d (debug mode)\n");
      printf("              -n iterations_thickness (10)\n");
      printf("              -i iterations_laplace (100)\n");
      printf("              -w (swap bytes of input) \n");
      printf("              -r (reverse) \n");
      printf("              -s (sum) \n");
      printf("              -l lambda (0.5)\n");
      printf("              -m compute and show mean thickness\n");
      printf("              -c color output (0: gray, 1: red-blue, 2: random)\n");
      printf("              --hx hy (1) \n");
      printf("              --hy hx (1) \n");
      printf("              --lw white matter label (3)\n");
      printf("              --lc cortex label (2)\n");
      printf("              --DT compute thicknes as DT\n");
      printf("              --streamlines show only the streamlines\n");
      return 1;
      break;
    
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  
  if ((argc - optind) != 2) {
    printf ("Incorrect number of arguments: ");
    printf("Author: Ruben Cardenes, April 2004 \n");
    printf("Usage: thickness2D [options] input2D.png output2D.png \n");
    printf("              (image dimensions are read from the input PNG)\n");
    printf("              -d (debug mode)\n");
    printf("              -n iterations_thickness (10)\n");
    printf("              -i iterations_laplace (100)\n");
    printf("              -w (swap bytes of input) \n");
    printf("              -r (reverse) \n");
    printf("              -s (sum) \n");
    printf("              -l lambda (0.5)\n");
    printf("              -m compute and show mean thickness\n");
    printf("              -c color output (0: gray, 1: red-blue, 2: random)\n");
    printf("              --hx hy (1) \n");
    printf("              --hy hx (1) \n");
    printf("              --lw white matter label (3)\n");
    printf("              --lc cortex label (2)\n");
    printf("              --DT compute thicknes as DT\n");
    printf("              --streamlines show only the streamlines\n");
    return 1;
  } else {
    if (sscanf(argv[optind++], "%s", inputfile) == 0)
      printf ("Error parsing argument \n");
    if (sscanf(argv[optind++], "%s", outputfile) == 0)
      printf ("Error parsing argument \n");
  }

  gettimeofday(&startinit,NULL);

  /* Read the domain from a PNG; the image dimensions become max1 (rows) and
     max2 (columns). */
  input = load_png_gray(inputfile, &width, &height);
  if (input == NULL) {
    exit(1);
  }
  max2 = width;
  max1 = height;
  printf("Input %s: %d rows x %d cols\n", inputfile, max1, max2);

  /* Report the distinct label values in the domain and require that the given
     --lw and --lc are actually present, so a wrong label fails loudly instead
     of silently producing an empty result. */
  {
    unsigned char present[256];
    print_domain_values(input, max1*max2, present);
    if (!require_label(present, label_wm, "--lw") ||
        !require_label(present, label_cortex, "--lc")) {
      free(input);
      exit(1);
    }
  }

  compute_boundary_cortex2D(input,max1,max2, label_wm, label_cortex);

  /* compute_boundary_cortex2D used the input --lw/--lc labels to detect the
     boundaries and then normalized the domain: the band is now labeled 2 (the
     interior 255, the inner boundary 1, the outer boundary 0). From here on the
     band is identified by that fixed label 2, regardless of the input --lc, so
     --lw/--lc may take any value present in the image. */
  label_cortex = 2;

  laplacefield  = (float**)malloc(sizeof(float*)*max1);
  for (i=0;i<max1;i++) {
    laplacefield[i] = (float*)malloc(sizeof(float)*max2);
  }

  printf("Entering in laplacian2D\n");
  if ( laplace2D(input, max1, max2, laplacefield, iterations_laplace, lambda, reverse) == 1 ) {
    printf("Error in thickness2D\n");
  }

  if (debug ==1) {
    printf("Writing ouput input_modificado.chr\n");
    if (debug == 1) {
      fp = fopen("input_modificado.chr","w"); 
      fwrite(input,sizeof(unsigned char),max1*max2,fp);     
      fclose(fp);
    }

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
  if (suma == 0) {
    if (reverse == 0) {
      if ( thickness2DYezzi(input, max1, max2, maps, laplacefield, 
			    gradientx, gradienty, num_it,hx,hy, label_cortex, debug) == 1 ) {
	printf("Error in thickness2D\n");
      }
    } else {
      if ( thickness2DYezzi_reverse(input, max1, max2, maps, laplacefield, 
				    gradientx, gradienty, num_it,hx,hy, label_cortex, debug) == 1 ) {
	printf("Error in thickness2D\n");
      }
    }
  }
  
  if (suma == 1) {
    maps_reverse = (float*)malloc(sizeof(float)*max1*max2);
    if ( thickness2DYezzi(input, max1, max2, maps, laplacefield, 
			  gradientx, gradienty, num_it,hx,hy, label_cortex, debug) == 1 ) {
      printf("Error in thickness2D\n");
    }
    if ( thickness2DYezzi_reverse(input, max1, max2, maps_reverse, laplacefield, 
				  gradientx, gradienty, num_it,hx,hy, label_cortex, debug) == 1 ) {
      printf("Error in thickness2D\n");
    }
    relabel_float(maps,max1*max2*max3,-1,0);
    relabel_float(maps_reverse,max1*max2*max3,-1,0);
    sumar_l1l2(maps,maps_reverse,maps,max1*max2*max3);     
    free(maps_reverse);
  }
  if (streamlines == 1) {
    streaml_maps = (float*)malloc(sizeof(float)*max1*max2);
    if ( thickness2Dgradient(input, max1, max2, maps, streaml_maps, laplacefield, gradientx, gradienty) == 1 ) {
      printf("Error in thickness2D\n");
    } 
    printf("Writing ouput streamline %s:\n",outputfile);
    write_float_output(outputfile, streaml_maps, max2, max1, 1, color_mode);
    free(streaml_maps);
  }


  gettimeofday(&endtotal,NULL);

  if (streamlines == 0 ) {
    printf("Writing ouput %s:\n",outputfile);
    write_float_output(outputfile, maps, max2, max1, 1, color_mode);
  }

  if (compute_mean == 1) {
    int npoints;
    float sigma, mean;
    mean = compute_mean_thickness2D(input, maps, label_cortex, max1, max2, &npoints, &sigma);
    printf("Mean thickness = %f std = %f (over %d band pixels)\n", mean, sigma, npoints);
  }

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
  return 0;
}
