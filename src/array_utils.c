#include "array_utils.h"

int relabel(unsigned char* data, int totdim, unsigned char oldlabel, unsigned char newlabel) {
  int i;
  for (i = 0; i < totdim; i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
}

int relabel_uchar(unsigned char* data, int totdim, unsigned char oldlabel, unsigned char newlabel) {
  int i;
  for (i = 0; i < totdim; i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
}

int relabel_ushort(unsigned short* data, int totdim, unsigned short oldlabel, unsigned short newlabel) {
  int i;
  for (i = 0; i < totdim; i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
}

int relabel_float(float* data, int totdim, float oldlabel, float newlabel) {
  int i;
  for (i = 0; i < totdim; i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
}

int sumar_l1l2(float* input1, float* input2, float* output, int totdim) {
  int i;

  for (i = 0; i < totdim; i++) {
    output[i] = input1[i] + input2[i] - 1;
  }

  return 0;
}
