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

  /* Report the distinct label values in the volume. laplace3D has no --lw/--lc
     to validate; it only needs a non-zero region over a zero background. */
  {
    unsigned char present[256];
    print_domain_values(input, max1*max2*max3, present);
  }

  output = (float***)malloc(sizeof(float**)*max3);
  for (k=0;k<max3;k++) {
    output[k] = (float**)malloc(sizeof(float*)*max1);
  }
  for (k=0;k<max3;k++) {
    for (i=0;i<max1;i++) {
      output[k][i] = (float*)malloc(sizeof(float)*max2);
    }
  }
  
  {
    unsigned char present[256];
    printf("Entering in EdgeDetect3D\n");
    if ( EdgeDetect3D(input, max1, max2, max3) == 1 ) {
      printf("Error in EdgeDetect3D\n");
      return 1;
    }
    printf("  after EdgeDetect3D: "); print_domain_values(input, max1*max2*max3, present);

    /* The Laplace field must be solved inside the domain material, so that region
       must carry label 2. After EdgeDetect3D the material interior keeps its
       original (non-zero) value and the surface is label 1; relabel every
       non-background, non-surface voxel to 2. RelabelBoundary3D then splits the
       two surfaces into the 0/1 Dirichlet values the solver reads. */
    printf("Relabeling\n");
    for (i=0;i<max1*max2*max3;i++) {
      if (input[i] != 0 && input[i] != 1) input[i] = 2;
    }
    printf("  after relabel(material->2): "); print_domain_values(input, max1*max2*max3, present);

    printf("Entering in RelabelBoundary3D\n");
    if ( RelabelBoundary3D(input, max1, max2, max3) == 1 ) {
      printf("Error in RelabelBoundary\n");
      return 1;
    }
    printf("  after RelabelBoundary3D: "); print_domain_values(input, max1*max2*max3, present);

    /* Label 2 is the only region the Laplace solver iterates over (see
       laplace3D()); every other label is held fixed. If the region you expect a
       smooth field in is not label 2, the solver never touches it. */
    {
      long counts[256]; int v;
      for (v=0;v<256;v++) counts[v]=0;
      for (i=0;i<max1*max2*max3;i++) counts[input[i]]++;
      printf("  label sizes going into solver (label:count):");
      for (v=0;v<256;v++) if (counts[v]) printf(" %d:%ld", v, counts[v]);
      printf("\n");
      printf("  --> solver iterates only over label 2: %ld voxels\n", counts[2]);
    }
  }
  
  printf("Writing domain domain_modificado3d.vol\n");
  fg=fopen("domain_modificado3d.vol","w");
  if (fg != NULL) {
    fwrite(input,sizeof(unsigned char),max1*max2*max3,fg);
    fclose(fg);
  } else {
    fprintf(stderr,"Failed writing domain_modificado3d.vol\n");
  }

  printf("Entering in laplacian3D\n");
  if ( laplace3D(input, max1, max2, max3, output, iterations, lambda) == 1 ) {
    printf("Error in laplace3D\n");

  }

  /* Output field statistics over the solved region (label 2). A healthy Laplace
     solution has a non-degenerate range there; a range of ~0, or 0 solved
     voxels, means the field is flat and nothing meaningful was computed. */
  {
    long nsolved=0; float smin=0,smax=0,ssum=0; int sinit=0, sum=0;
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
        for (j=0;j<max2;j++) {
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

