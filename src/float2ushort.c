#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

long fileSize(char *fname)
{
  struct stat statBuf;

  if (stat(fname, &statBuf) == -1) {
    /* fprintf(stderr,"Failed to stat %s\n",fname); */
    return -1;
  }
  return statBuf.st_size;
}

int main(int argc,char** argv){
  FILE *fp = NULL;
  float *data = NULL;
  unsigned short *dataout;
  int npixels = -1;
  int i;

  if (argc != 3) { 
    printf("Usage: float2ushort filein fileout\n");
    return; 
  } 
  fp = fopen(argv[1],"r");
  if (fp == NULL) {
    printf("No such file \n");
    return;
  }
  npixels = (int)fileSize(argv[1])/4;
  printf("size: %d \n", npixels);
  data = (float *)malloc(sizeof(float)*npixels);
  if (fread(data, sizeof(float), npixels, fp) == 0) {
    fprintf(stderr,"Failed reading the data\n");
    fclose(fp);
    free(data);
    return -1;
  }

  dataout = (unsigned short *)malloc(sizeof(unsigned short)*npixels);
  for (i=0;i<npixels;i++) {
    dataout[i]=(unsigned short)data[i];
  }
  
  fclose(fp);
  fp = fopen(argv[2],"w");
  fwrite(dataout,sizeof(unsigned short),npixels,fp);

  free(data);
  free(dataout);

  printf("float2ushort OK\n");
}









