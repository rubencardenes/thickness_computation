/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>   /* getopt, optarg, optind (with a proper prototype) */
#include "laplace2D.h"
#include "png_write.h"

int main(int argc,char* argv[]) {
  unsigned char* input;
  float **output;
  int i;
  int iterations = 10;
  float lambda = 0.5;
  int color_mode = COLOR_GRAY;
  int width, height, c, nargs;
  char *inputfile,*outputfile;
  FILE *fg;

  while ((c = getopt(argc, argv, "c:")) != -1) {
    switch (c) {
    case 'c': color_mode = atoi(optarg); break;
    default: break;
    }
  }

  nargs = argc - optind;
  if (nargs != 3 && nargs != 4) {
    printf("Usage: laplace2D [-c color] input.png output.png iterations [lambda]\n");
    printf("       (image dimensions are read from the input PNG)\n");
    printf("       -c color output (0: gray, 1: red-blue, 2: random)\n");
    return 1;
  }

  inputfile  = argv[optind + 0];
  outputfile = argv[optind + 1];
  iterations = atoi(argv[optind + 2]);
  if (nargs == 4) {
    lambda = atof(argv[optind + 3]);
  }

  /* Read the domain from a PNG; dimensions become height (rows) and width (cols). */
  input = load_png_gray(inputfile, &width, &height);
  if (input == NULL) {
    exit(1);
  }
  printf("Input %s: %d rows x %d cols\n", inputfile, height, width);

  /* Report the distinct label values in the domain. laplace2D has no --lw/--lc
     to validate; it only needs a non-zero region over a zero background. */
  {
    unsigned char present[256];
    print_domain_values(input, height*width, present);
  }

  output = (float**)malloc(sizeof(float*)*height);
  for (i=0;i<height;i++) {
    output[i] = (float*)malloc(sizeof(float)*width);
  }
  
  {
    unsigned char present[256];
    printf("Entering in EdgeDetect\n");
    if ( EdgeDetect(input, height, width) == 1 ) {
      printf("Error in EdgeDetect\n");
    }
    printf("  after EdgeDetect: "); print_domain_values(input, height*width, present);

    /* The Laplace field must be solved inside the domain material (the ring
       wall), so that region must carry label 2. After EdgeDetect the material
       interior keeps its original (non-zero) value and the edges are label 1;
       relabel every non-background, non-edge pixel to 2. RelabelBoundary then
       splits the two edges into the 0/1 Dirichlet values the solver reads. */
    printf("Relabeling\n");
    for (i=0;i<height*width;i++) {
      if (input[i] != 0 && input[i] != 1) input[i] = 2;
    }
    printf("  after relabel(material->2): "); print_domain_values(input, height*width, present);

    printf("Entering in RelabelBoundary\n");
    if ( RelabelBoundary(input, height, width) == 1 ) {
      printf("Error in RelabelBoundary\n");
    }
    printf("  after RelabelBoundary: "); print_domain_values(input, height*width, present);

    /* Report the size of each label class going into the solver. Label 2 is the
       only region the Laplace solver iterates over (see laplace2D()); every other
       label is held fixed at its value. If the region you expect to see a smooth
       field in is not label 2, the solver never touches it. */
    {
      long counts[256]; int v;
      for (v=0;v<256;v++) counts[v]=0;
      for (i=0;i<height*width;i++) counts[input[i]]++;
      printf("  label sizes going into solver (label:count):");
      for (v=0;v<256;v++) if (counts[v]) printf(" %d:%ld", v, counts[v]);
      printf("\n");
      printf("  --> solver iterates only over label 2: %ld pixels\n", counts[2]);
    }
  }
  
  printf("Writing domain domain_anillo_modificado.chr\n");
  fg=fopen("domain_anillo_modificado.chr","w");
  if (fg != NULL) {
    fwrite(input,sizeof(unsigned char),height*width,fg);
    fclose(fg);
  } else {
    fprintf(stderr,"Failed writing domain_anillo_modificado.chr\n");
  }

  printf("Entering in laplacian2D\n");
  if ( laplace2D(input, height, width, output, iterations, lambda, 0) == 1 ) {
    printf("Error in thickness2D\n");
  }

  /* Output field statistics, split between the solved region (label 2) and the
     rest. A healthy Laplace solution has a non-degenerate range over the solved
     region (values spread between the two boundary conditions). A range of ~0
     there, or 0 solved pixels, means the field is flat and nothing meaningful
     was computed. */
  {
    int r,col; long nsolved=0, nother=0;
    float smin=0,smax=0,ssum=0, omin=0,omax=0;
    int sinit=0, oinit=0;
    for (r=0;r<height;r++) {
      for (col=0;col<width;col++) {
        float v = output[r][col];
        if (input[r*width+col] == 2) {
          if (!sinit) { smin=smax=v; sinit=1; } else { if (v<smin) smin=v; if (v>smax) smax=v; }
          ssum += v; nsolved++;
        } else {
          if (!oinit) { omin=omax=v; oinit=1; } else { if (v<omin) omin=v; if (v>omax) omax=v; }
          nother++;
        }
      }
    }
    printf("Field stats: solved region (label 2): %ld px, min=%.4f max=%.4f mean=%.4f (range=%.4f)\n",
           nsolved, smin, smax, nsolved?ssum/nsolved:0.0f, smax-smin);
    printf("             everything else: %ld px, min=%.4f max=%.4f (held fixed at input labels)\n",
           nother, omin, omax);
    if (nsolved == 0) {
      printf("  WARNING: no pixels labelled 2, so the Laplace solver did nothing.\n");
    } else if (smax-smin < 1e-6f) {
      printf("  WARNING: solved field is flat (range ~0); check boundary labels.\n");
    }
  }

  printf("Writing ouput %s:\n",outputfile);
  {
    /* Copy the row-pointer field into a contiguous buffer so it can be written
       as a raw float file or a min..max normalized (optionally colored) PNG. */
    float *buf = (float*)malloc(sizeof(float)*height*width);
    int r,col;
    for (r=0;r<height;r++) {
      for (col=0;col<width;col++) {
        buf[r*width+col] = output[r][col];
      }
    }
    /* For the PNG, show only the solved Laplacian field (label 2). The other
       labels are held at fixed values (0/1 boundaries, 3 background) that would
       otherwise dominate the min..max normalization and crush the field's range;
       masking them to NaN makes them render black and excludes them from the
       normalization, so the field itself is stretched across the full 0..255.
       The raw .flt output keeps the complete field unchanged. */
    if (png_has_extension(outputfile)) {
      for (r=0;r<height*width;r++) {
        if (input[r] != 2) buf[r] = NAN;
      }
    }
    write_float_output(outputfile, buf, width, height, 0, color_mode);
    free(buf);
  }

  free(input);
  free(output);
  printf("OK laplace2D\n");
  return 0;
}

