/* Copyright (c) Ruben Cardenes Almeida 15/03/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include "io.h"
#include "compute_boundary_cortex3D.h"

int main(int argc, char *argv[]) {
  unsigned short *input;
  int i, swapbyte = 0, label_wm = 3, label_cortex = 2, c, option_index;
  int height = 256;
  int width = 256;
  int depth = 1;
  int debug = 1;
  FILE *fp;
  struct timeval startinit;
  char input_prefix[200];
  char input_file[200], outputfile[200], extension[6];
  unsigned char *aux, *l;

  while (1) {
    static struct option long_options[] = {
        {0, 0, 0, 0}};

    c = getopt_long(argc, argv, "dwl:c:", long_options, &option_index);

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
      printf("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if ((argc - optind) != 5) {
    printf("Incorrect number of arguments: ");
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
        printf("Error parsing argument \n");
      if (sscanf(argv[optind++], "%s", outputfile) == 0)
        printf("Error parsing argument \n");
      if (sscanf(argv[optind++], "%d", &height) == 0)
        printf("Error parsing argument \n");
      if (sscanf(argv[optind++], "%d", &width) == 0)
        printf("Error parsing argument \n");
      if (sscanf(argv[optind++], "%d", &depth) == 0)
        printf("Error parsing argument \n");
    }
  }

  gettimeofday(&startinit, NULL);
  /* reserve data */
  input = (unsigned short *)malloc(sizeof(unsigned short) * height * width * depth);

  l = strchr(input_prefix, '.');
  if (l != (unsigned char *)NULL) {
    strncpy(extension, (char *)l, 5);
    extension[5] = '\0';
  } else {
    strcpy(extension, "");
  }
  if (strcmp(extension, ".vols") == 0) {
    fp = fopen(input_prefix, "r");
    if (fp == NULL) {
      fprintf(stderr, "Failed reading inputfile %s\n", input_prefix);
      exit(1);
    }
    fread(input, sizeof(unsigned short), height * width * depth, fp);
    fclose(fp);
  } else {
    /* Read data */
    for (i = 0; i < depth; i++) {
      sprintf(input_file, "%s.%03d", input_prefix, i + 1);
      fp = fopen(input_file, "r");
      if (fp == NULL) {
        fprintf(stderr, "Failed reading inputfile %s\n", input_file);
        exit(1);
      }
      fread(&input[height * width * i], sizeof(unsigned short), height * width, fp);
      fclose(fp);
    }
  }

  if (swapbyte == 1) {
    for (i = 0; i < height * width * depth; i++) {
      aux = (unsigned char *)&input[i];
      input[i] = ReadGEShort(aux);
    }
  }

  /* Compute */
  compute_boundary_cortex3D(input, height, width, depth, label_wm, label_cortex);

  /* Write data */
  printf("writing file %s\n", outputfile);
  fp = fopen(outputfile, "w");
  if (fp == NULL) {
    fprintf(stderr, "Failed writing file %s\n", outputfile);
    return 1;
  }
  fwrite(input, sizeof(unsigned short), height * width * depth, fp);
  fclose(fp);

  /* Free data */
  free(input);
  printf("compute_boundary_cortex3D OK\n");
  return 0;
}
