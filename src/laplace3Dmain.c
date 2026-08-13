/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "laplace3D.h"
#include "png_write.h"

int main(int argc,char* argv[]) {
  unsigned char* input;
  float ***output;
  int i,j,k;
  int height = 256;
  int width = 256;
  int depth = 70;
  int iterations = 10;
  float lambda = 0.5;
  FILE *fp,*fg;
  if (argc != 8 && argc != 7) {
    printf("Usage: laplace3D height width depth input.vol output.volf iterations [lambda]\n");
    return 1;
  }

  height = atoi(argv[1]);
  width = atoi(argv[2]);
  depth = atoi(argv[3]);
  iterations = atoi(argv[6]);
  if (argc == 8) {
    lambda = atof(argv[7]);
  }
  
  input = (unsigned char*)malloc(sizeof(unsigned char)*height*width*depth);
  fp = fopen(argv[4],"rb");
  if (fp == NULL) {
    fprintf(stderr,"Failed reading inputfile %s\n",argv[4]);
    exit(1);
  }
  fread(input,sizeof(unsigned char),height*width*depth,fp);
  fclose(fp);

  /* Report the distinct label values in the volume. laplace3D has no --lw/--lc
     to validate; it only needs a non-zero region over a zero background. */
  {
    unsigned char present[256];
    print_domain_values(input, height*width*depth, present);
  }

  output = (float***)malloc(sizeof(float**)*depth);
  for (k=0;k<depth;k++) {
    output[k] = (float**)malloc(sizeof(float*)*height);
  }
  for (k=0;k<depth;k++) {
    for (i=0;i<height;i++) {
      output[k][i] = (float*)malloc(sizeof(float)*width);
    }
  }
  
  {
    unsigned char present[256];
    printf("Entering in EdgeDetect3D\n");
    if ( EdgeDetect3D(input, height, width, depth) == 1 ) {
      printf("Error in EdgeDetect3D\n");
      return 1;
    }
    printf("  after EdgeDetect3D: "); print_domain_values(input, height*width*depth, present);

    /* The Laplace field must be solved inside the domain material, so that region
       must carry label 2. After EdgeDetect3D the material interior keeps its
       original (non-zero) value and the surface is label 1; relabel every
       non-background, non-surface voxel to 2. RelabelBoundary3D then splits the
       two surfaces into the 0/1 Dirichlet values the solver reads. */
    printf("Relabeling\n");
    for (i=0;i<height*width*depth;i++) {
      if (input[i] != 0 && input[i] != 1) input[i] = 2;
    }
    printf("  after relabel(material->2): "); print_domain_values(input, height*width*depth, present);

    printf("Entering in RelabelBoundary3D\n");
    if ( RelabelBoundary3D(input, height, width, depth) == 1 ) {
      printf("Error in RelabelBoundary\n");
      return 1;
    }
    printf("  after RelabelBoundary3D: "); print_domain_values(input, height*width*depth, present);

    /* Label 2 is the only region the Laplace solver iterates over (see
       laplace3D()); every other label is held fixed. If the region you expect a
       smooth field in is not label 2, the solver never touches it. */
    {
      long counts[256]; int v;
      for (v=0;v<256;v++) counts[v]=0;
      for (i=0;i<height*width*depth;i++) counts[input[i]]++;
      printf("  label sizes going into solver (label:count):");
      for (v=0;v<256;v++) if (counts[v]) printf(" %d:%ld", v, counts[v]);
      printf("\n");
      printf("  --> solver iterates only over label 2: %ld voxels\n", counts[2]);
    }
  }
  
  printf("Writing domain domain_modificado3d.vol\n");
  fg=fopen("domain_modificado3d.vol","w");
  if (fg != NULL) {
    fwrite(input,sizeof(unsigned char),height*width*depth,fg);
    fclose(fg);
  } else {
    fprintf(stderr,"Failed writing domain_modificado3d.vol\n");
  }

  printf("Entering in laplacian3D\n");
  if ( laplace3D(input, height, width, depth, output, iterations, lambda) == 1 ) {
    printf("Error in laplace3D\n");

  }

  /* Output field statistics over the solved region (label 2). A healthy Laplace
     solution has a non-degenerate range there; a range of ~0, or 0 solved
     voxels, means the field is flat and nothing meaningful was computed. */
  {
    long nsolved=0; float smin=0,smax=0,ssum=0; int sinit=0, sum=0;
    for (k=0;k<depth;k++) {
      for (j=0;j<width;j++) {
        for (i=0;i<height;i++) {
          if (input[sum] == 2) {
            float v = output[k][i][j];
            if (!sinit) { smin=smax=v; sinit=1; } else { if (v<smin) smin=v; if (v>smax) smax=v; }
            ssum += v; nsolved++;
          }
          sum++;
        }
      }
    }
    printf("Field stats: solved region (label 2): %ld vox, min=%.4f max=%.4f mean=%.4f (range=%.4f)\n",
           nsolved, smin, smax, nsolved?ssum/nsolved:0.0f, smax-smin);
    if (nsolved == 0) {
      printf("  WARNING: no voxels labelled 2, so the Laplace solver did nothing.\n");
    } else if (smax-smin < 1e-6f) {
      printf("  WARNING: solved field is flat (range ~0); check boundary labels.\n");
    }
  }

  printf("Writing ouput %s:\n",argv[5]);
  fg=fopen(argv[5],"w");
  for (k=0;k<depth;k++) {
    for (i=0;i<height;i++) {
      fwrite(output[k][i],sizeof(float),width,fg);    
    }
  }
  fclose(fg);

  free(input);
  free(output);
  printf("OK laplace2D\n");
  return 0;
}

