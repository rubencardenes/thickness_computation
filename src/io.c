#include <sys/stat.h>
#include "io.h"

unsigned short ReadGEShort(unsigned char *x)
// Convert from GE short format to short on machine
// independent of byte order
{
  return (x[0] * (1 << 8) + x[1]);
  // MSB               LSB  in GE format
}

/* FILE OPERATIONS */
long fileSize(char *fname) {
  struct stat statBuf;

  if (stat(fname, &statBuf) == -1) {
    /* fprintf(stderr,"Failed to stat %s\n",fname); */
    return -1;
  }
  return statBuf.st_size;
}

/* This basically assumes I only want to read 1024^2, 512^2 or 256^2 images */
long headerSize(int fileSize) {
  long headerSize = 0;

  /* Assume it is a 512^2 genesis format image */
  if (fileSize == 532192) {
    headerSize = 7904;
  } else if (fileSize == 138976) {
    headerSize = 7904;
  } else if (fileSize == 145408) {
    /* 256x256 SIGNA header */
    headerSize = 14336;
  } else if (fileSize >= 1024 * 1024 * 2) {
    headerSize = fileSize - 1024 * 1024 * 2;
  } else if (fileSize >= 768 * 768 * 2) {
    headerSize = fileSize - 768 * 768 * 2;
  } else if (fileSize >= 512 * 512 * 2) {
    headerSize = fileSize - 512 * 512 * 2;
  } else if (fileSize >= 320 * 320 * 2) {
    headerSize = fileSize - 320 * 320 * 2;
  } else if (fileSize >= 256 * 256 * 2) {
    headerSize = fileSize - 256 * 256 * 2;
  }
  return headerSize;
}
