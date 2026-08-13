#ifndef IO_H
#define IO_H

#include <stdio.h>

unsigned short ReadGEShort(unsigned char *x);
long fileSize(char *fname);
long headerSize(int fileSize);
#endif /* IO_H */
