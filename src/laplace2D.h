
/* Gradient definitions */
#define FDY(Data,iX,iY) (Data[(iX)][(iY)+1] - Data[(iX)][(iY)])
#define BDY(Data,iX,iY) (Data[(iX)][(iY)] - Data[(iX)][(iY)-1])

#define FDX(Data,iX,iY) (Data[(iX)+1][(iY)] - Data[(iX)][(iY)])
#define BDX(Data,iX,iY) (Data[(iX)][(iY)] - Data[(iX)-1][(iY)])

#define CDX(Data,iX,iY) ((Data[(iX)+1][(iY)] - Data[(iX)-1][(iY)])*0.5)
#define CDY(Data,iX,iY) ((Data[(iX)][(iY)+1] - Data[(iX)][(iY)-1])*0.5)


int maptox(int mapindex,int max2);
int maptoy(int mapindex,int max2);
float distance(int x1,int y1,int x2,int y2);
int compute_corners(unsigned short *input, int max1,int max2);
int new_compute_corners(unsigned short *input, int max1,int max2);
int relabel(unsigned char* data, int totdim, unsigned char oldlabel, unsigned char newlabel);
int relabel_ushort(unsigned short* data, int totdim, int oldlabel, int newlabel);
int relabel_float(float* data, int totdim, int oldlabel, int newlabel);
int laplace2D(unsigned char* input,int max1, int max2, float** output, int iterations, float lambda, int reverse);
int RelabelBoundary(unsigned char *domain,int max1,int max2);
int EdgeDetect(unsigned char *domain, int max1, int max2);
int normalize(float** gradientx, float** gradienty, int max1,int max2);
int iGradX(float **ppfData, float **ppfGradient, int numRowX, int numColY);
int iGradY(float **ppfData, float **ppfGradient, int numRowX, int numColY);
int maxcomponent2D(unsigned short* data, int max1,int max2, int label);
int sizefilter2D(unsigned short* data, int max1,int max2, int tam, int oldlabel, int newlabel);
int sumar_l1l2(float* input1, float* input2, float* output, int totdim);
