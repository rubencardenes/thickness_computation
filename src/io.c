#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include "io.h"

int ReadFromFile(void *ptr, int size, int nitems,int headerSize, FILE *fp) {
  if (headerSize > 0) {
    if (fseek(fp, headerSize, 1) != 0) return -1;
  }
  return fread(ptr, size, nitems,fp);
}

int IsMSBFirstForShort(void)
/* Return 1 if MSB of a short is at the first byte */
{
  union {
    short i;
    char c[sizeof(short)];
  } kludge;
  kludge.i = 0;
  kludge.c[0] = 1;
 
  if (kludge.i != 1) {
    return 1;
  } else {
    return 0;
  }
}

unsigned short ReadGEShort(unsigned char *x)
// Convert from GE short format to short on machine
// independent of byte order
{
  return (x[0] * (1 << 8) + x[1]);
      // MSB               LSB  in GE format
}

// All functions to read different file formats
int ReadRawMRI(unsigned short *ptr, int headerSize, FILE *fp, int Nrow, int Ncol)
{
  int i,j;
  unsigned char* aux;

  if (headerSize > 0) {
    if (fseek(fp, headerSize, 1) != 0) return -1;
  }

  int iBytesRead;
  iBytesRead = fread(ptr, sizeof(unsigned short), Nrow*Ncol, fp);
  if ((iBytesRead != Nrow*Ncol) || ferror(fp) ) {
    printf("Error: Failed while reading data\n");
    free(ptr);
    return -1;
  }

  // The GE short format has MSB first and LSB second
  // Some machines (i486) use LSB first and MSB second
  // so if this machine doesn't have the MSB first we must convert it
  if (! IsMSBFirstForShort()) {
    // Convert each short in the array
    for (i=0; i < Nrow*Ncol; i++) {
      aux = (unsigned char*)&ptr[i];
      ptr[i] = ReadGEShort(aux);
    }
  }
  return 0;
}

int WriteMRI(unsigned short *data,int headersize,int nrow, int ncol, FILE* fp,char* headerFName) {
  FILE *header = (FILE *)NULL;
  unsigned char *byteHeader = 0;

  if (headersize == 0) {
    fwrite(data,sizeof(unsigned short),nrow*ncol,fp);
    return 0;
  } 
  
  if (headersize != 0) {
    if ((header = fopen(headerFName,"rb")) == NULL) return 6;
  
    byteHeader = (unsigned char*)malloc(sizeof(unsigned char)*headersize);
    if (byteHeader == (unsigned char *)NULL) {
      return 1; // Out of memory error
    }
    // Read the header from the file header
    if (fread(byteHeader,sizeof(unsigned char),headersize,header) != headersize) {
      // fread returns 0 on error or eof
      free(byteHeader);
      return 1;
    }
    fclose(header);
  
    if (fwrite(byteHeader,sizeof(unsigned char),headersize,fp) != headersize) {
      printf("Error during writing\n");
    }

    free(byteHeader);
    fwrite(data,sizeof(unsigned short),nrow*ncol,fp);
    return 0;
  }

}

/* FILE OPERATIONS */
long fileSize(char *fname)
{
  struct stat statBuf;

  if (stat(fname, &statBuf) == -1) {
    /* fprintf(stderr,"Failed to stat %s\n",fname); */
    return -1;
  }
  return statBuf.st_size;
}

/* This basically assumes I only want to read 1024^2, 512^2 or 256^2 images */
long headerSize(int fileSize)
{
  long headerSize = 0;

  /* Assume it is a 512^2 genesis format image */
  if (fileSize == 532192) {
    headerSize = 7904;
  } else if (fileSize == 138976) {
    headerSize = 7904;
  } else if (fileSize == 145408) {
    /* 256x256 SIGNA header */
    headerSize = 14336;
  } else if (fileSize >= 1024*1024*2) {
    headerSize = fileSize - 1024*1024*2;
  } else if (fileSize >= 768*768*2) {
    headerSize = fileSize - 768*768*2;
  } else if (fileSize >= 512*512*2) {
    headerSize = fileSize - 512*512*2;
  } else if (fileSize >= 320*320*2) {
    headerSize = fileSize - 320*320*2;
  } else if (fileSize >= 256*256*2) {
    headerSize = fileSize - 256*256*2;
  }
  return headerSize;
}
