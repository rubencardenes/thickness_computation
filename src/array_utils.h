#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

int relabel(unsigned char* data, int totdim, unsigned char oldlabel, unsigned char newlabel);
int relabel_ushort(unsigned short* data, int totdim, unsigned short oldlabel, unsigned short newlabel);
int relabel_float(float* data, int totdim, float oldlabel, float newlabel);
int sumar_l1l2(float* input1, float* input2, float* output, int totdim);

#endif /* ARRAY_UTILS_H */
