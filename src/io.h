#include <stdio.h>

int ReadFromFile(void *ptr, int size, int nitems, int headerSize, FILE *fp);
int IsMSBFirstForShort(void);
unsigned short ReadGEShort(unsigned char *x);
int ReadRawMRI(unsigned short *ptr, int headerSize, FILE *fp, int Nrow, int Ncol);
int WriteMRI(unsigned short *data,int headersize,int nrow, int ncol, FILE* fp,char* headerFName);
long fileSize(char *fname);
long headerSize(int fileSize);
