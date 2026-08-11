/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "laplace2D.h"

typedef  unsigned char uchar;
#define SEVEN_SUPP     
#define MAX_CORNERS   15000  /* max corners per frame */
#define  FTOI(a) ( (a) < 0 ? ((int)(a-0.5)) : ((int)(a+0.5)) )
typedef  struct {int x,y,info, dx, dy, I;} CORNER_LIST[MAX_CORNERS];

struct fifo{
  int *index_elem;
  int first;
  int last;
};

int maptox(int mapindex,int max2) {
  return mapindex % max2;
}

int maptoy(int mapindex,int max2) {
  return (mapindex - (mapindex % max2)) / max2;
}

float distance(int x1,int y1,int x2,int y2) {
  return sqrt ((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
}

/* Input values
255 Outside domain
1 Exterior boundary 
0 Interior boundary
2 Inside domain 
*/ 
int laplace2D(unsigned char* input,int max1, int max2, float** output, int iterations, float lambda, int reverse) {
  int i,j,l;
  int sum = 0;
  /* Initialize domain, inside=0, and boundaries values*/
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (input[sum] == 2) {
	output[i][j] = 0;
      } else {
	if (input[sum] == 255) {
	  if (reverse) {
	    output[i][j] = -1;
	  } else {
	    output[i][j] = 255;
	  }
	} else {
	  output[i][j] = input[sum];
	}
      }
      sum++;
    }
  }
  
  /* Solve Laplacian */
  for (l=0;l<iterations;l++) {
    sum = 0;
    for (i=0;i<max1;i++) {
      for (j=0;j<max2;j++) {
	if (input[sum] == 2 
	    && i != 0 && i != max1-1 && j!=0 && j != max2-1) {
	  /*output[i][j] = 0.25 *(output[i-1][j] + output[i+1][j] + output[i][j-1] + output[i][j+1]);*/
	  output[i][j] = output[i][j]+(lambda+1)*(0.25 *(output[i-1][j] + output[i+1][j] + output[i][j-1] + output[i][j+1]) - output[i][j]);
	}
	sum++;
      }
    }
  }
 
  return 0;
}

int EdgeDetect(unsigned char *domain, int max1, int max2){
  int x,y,i;
  i = 0;
 
  for(x=0; x<max1; x++) {
    for(y=0; y<max2; y++) {
      if ((x==0)||(y==0)||(x==max1-1)||(y==max2-1)) {
	/* domain[i]=255;*/
      } else
	if ( (domain[i]!=0)&&
	     ((domain[i+1]==0)||
	      (domain[i-1]==0)||

	      (domain[i+max2]==0)||
	      (domain[i-max2]==0))) {

	  domain[i]=1;
	}
      /*else {	   	   
	domain[i]=255;
	}*/
      i++;
    }
  }
  return 0;
}


int maxcomponent2D(unsigned short* data, int max1,int max2, int label) {
  int i,j,k,xr,yr,x,y,mapindex,newmapindex,exceed,count,xnew,ynew;
  struct list {
    int num_elem;
    int *elem;
  };
  struct list mylist,list2;
  int max_num_elem = max1*max2;
  unsigned char *aux_data; 
  int tam = 0;

  aux_data = (unsigned char*)calloc(max1*max2,sizeof(unsigned char));
  mylist.elem = (int*)malloc(sizeof(int)*max_num_elem);
  mylist.num_elem = 0;
  list2.elem = (int*)malloc(sizeof(int)*max_num_elem);
  list2.num_elem = 0; 

  count = 0;
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (data[count] == label && aux_data[count] == 0) {
	mylist.elem[mylist.num_elem] = i*max2+j;
	mylist.num_elem++;
	list2.elem[list2.num_elem] = i*max2+j;
	list2.num_elem++;    
	while (mylist.num_elem !=0) {
	  mapindex = mylist.elem[mylist.num_elem-1];
	  mylist.num_elem--;
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      newmapindex = mapindex + max2*y + x;
	      xnew = maptox(newmapindex,max2);
	      ynew = maptoy(newmapindex,max2);
	      if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 &&
		  data[newmapindex] == label && aux_data[newmapindex] == 0) {
		mylist.elem[mylist.num_elem] = newmapindex;
		mylist.num_elem++;
		list2.elem[list2.num_elem] = newmapindex;
		list2.num_elem++;
		aux_data[newmapindex] = 1;
	      }
	    }
	  }
	}
	if (list2.num_elem > tam) {
	  tam = list2.num_elem;
	}
	list2.num_elem = 0;
	mylist.num_elem = 0;
      }
      count++;
    }
  }

  free(mylist.elem);
  free(list2.elem);
  free(aux_data);

  return tam;

}

int sizefilter2D(unsigned short* data, int max1,int max2, int tam, int oldlabel, int newlabel) {
  int i,j,k,xr,yr,x,y,mapindex,newmapindex,exceed,count,xnew,ynew;
  struct list {
    int num_elem;
    int *elem;
  };
  struct list mylist,list2;
  int max_num_elem = max1*max2;
  unsigned char *aux_data; 

  aux_data = (unsigned char*)calloc(max1*max2,sizeof(unsigned char));
  mylist.elem = (int*)malloc(sizeof(int)*max_num_elem);
  mylist.num_elem = 0;
  list2.elem = (int*)malloc(sizeof(int)*max_num_elem);
  list2.num_elem = 0; 

  count = 0;
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (data[count] == oldlabel && aux_data[count] == 0) {
	mylist.elem[mylist.num_elem] = i*max2+j;
	mylist.num_elem++;
	list2.elem[list2.num_elem] = i*max2+j;
	list2.num_elem++;    
	exceed = 0;
	while (mylist.num_elem !=0) {
	  mapindex = mylist.elem[mylist.num_elem-1];
	  mylist.num_elem--;
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      newmapindex = mapindex + max2*y + x;
	      xnew = maptox(newmapindex,max2);
	      ynew = maptoy(newmapindex,max2);
	      if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 &&
		  data[newmapindex] == oldlabel && aux_data[newmapindex] == 0) {
		mylist.elem[mylist.num_elem] = newmapindex;
		mylist.num_elem++;
		list2.elem[list2.num_elem] = newmapindex;
		list2.num_elem++;
		aux_data[newmapindex] = 1;
		if (list2.num_elem > tam) {
		  exceed = 1;
		}
	      }
	    }
	  }
	}
	if (exceed == 0) { /* If the size is not exceeded we relabel de data */	 
	  for (k=0;k<list2.num_elem;k++) {
	    mapindex = list2.elem[k];
	    data[mapindex] = newlabel;
	  }  
	}
	list2.num_elem = 0;
	mylist.num_elem = 0;
      }
      count++;
    }
  }

  free(mylist.elem);
  free(list2.elem);
  free(aux_data);

  return 0;

}

int relabel_float(float* data, int totdim, int oldlabel, int newlabel) {
  int i;

  for (i=0;i<totdim;i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }

  return 0;
}

int relabel_ushort(unsigned short* data, int totdim, int oldlabel, int newlabel) {
  int i;

  for (i=0;i<totdim;i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }

  return 0;
}

int relabel(unsigned char* data, int totdim, unsigned char oldlabel, unsigned char newlabel) {
  int i;
  for (i=0;i<totdim;i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
} 

int floodfill(unsigned char *domain, int startindex, unsigned short oldlabel, unsigned char newlabel,int max1, int max2) {
  int x,y,newmapindex,xnew,ynew;
  if (domain[startindex] == oldlabel) {
    domain[startindex] = newlabel;
  }
  for (x=-1;x<2;x++) {
    for (y=-1;y<2;y++) {
      newmapindex = startindex + max2*y + x;
      xnew = maptox(newmapindex,max2);
      ynew = maptoy(newmapindex,max2);
      if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1) {
	if (domain[newmapindex] == oldlabel) {
	  floodfill(domain,newmapindex,oldlabel,newlabel,max1,max2);
	}
      }
    }
  }
}

int RelabelBoundary(unsigned char *domain,int max1,int max2){
  int x,y,xr,yr,i,xnew,ynew;
  int newmapindex,mapindex,start,newstart;
  struct list {
    int num_elem;
    int *elem;
  } mylist;
  int max_num_elem_mylist = max1*max2;

  i = 0;
  start = -1;

  for(y=0; y<max2; y++) {
    for(x=0; x<max1; x++) {
      if (domain[i] == 1) {
	start = i;
	y=max1;x=max2;
      }
      i++;
    }
  }

  /* No boundary pixel (value 1) was found: the input domain does not have the
     expected form (a non-zero region bordering a zero background), so there is
     nothing to relabel. Bail out instead of writing to an invalid index. */
  if (start < 0) {
    fprintf(stderr,"RelabelBoundary: no boundary pixel found; check input labels\n");
    return 1;
  }

  mylist.elem = (int*)malloc(sizeof(int)*max_num_elem_mylist);
  domain[start] = 0;
  mapindex = start;
  mylist.num_elem = 1;
  mylist.elem[0] = start;
  while (mylist.num_elem != 0) {
    /* Get new element from mylist */
    mapindex = mylist.elem[mylist.num_elem-1];
    mylist.num_elem--;
    for (x=-1;x<2;x++) {
      for (y=-1;y<2;y++) {
	if (x==0 && y ==0) continue;
	newmapindex = mapindex + x + y*max2;
	xnew = maptox(newmapindex,max2);
	ynew = maptoy(newmapindex,max2);
	if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 && domain[newmapindex] == 1) {
	  /* Put new element in mylist*/
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  if (mylist.num_elem > max_num_elem_mylist) {
	    printf("Error, exceeded number of elements in mylist.num_elem ( > %d ), in RelabelBoundary\n",max_num_elem_mylist);
	    return 1;
	  }
	  domain[newmapindex] = 0;
	}
      }
    }
  }

  /* find the seed for the floodfill */
  /* for(y=0; y<max2; y++) {
    for(x=0; x<max1; x++) {
      if (domain[i] == 1) {
	start = i;
	for (xr=-1;xr<2;xr++) {
	  for (yr=-1;yr<2;yr++) {
	    if (domain[start+yr*max2+xr] == 255) {
	      newstart = start+yr*max2+xr;
	      yr = 2; xr =2;
	    }	  
	  }
	}
	start = i;
	y=max1;x=max2;
      }
      i++;    
    }
    }*/
  
  /* floodill the interior region, relabeling to -1 */ 
  floodfill(domain,0,domain[0],3,max1,max2);

  free(mylist.elem);
  return 0;
}


int iGradY(float **ppfData, float **ppfGradient, int numRowX, int numColY)
/* 
 * PURPOSE : Calculate the partial derivative in X direction.
 * 			X - increases across the rows.
 * 			Y - increases across the columns.
 * INPUTS :
 * ppfData     : The image data.
 * ppfGradient : The gradient will be stored in this array
 * numRowX	   : The number of points in the X direction.
 * numColY	   : The number of points in the Y direction.
 *
 * OUTPUTS :
 * 		Ix is approximated.
 *
 * METHOD : At each point in the image a central difference approx.
 * 	    to the first derivative in the x direction is made.
 * 	    At the boundaries a forward difference or a backward difference
 * 	    approximation is made.
 */
{
    register int i; register int j; register int iMaxX; register int iMaxY;

    iMaxX = numRowX - 1;
    iMaxY = numColY - 1;
    assert(iMaxX >= 0 && iMaxY >= 0);

    if (iMaxX == 0) { /* Ix == 0 everywhere for one row data */
	  for (j = 0; j <= iMaxY; j++) {
		ppfGradient[0][j] = 0;
	  }
    } else { /* iMaxY can still be zero - but this is OK */
		// Along edges of image where i == 0 and i == iMaxX (top and bottom row)
		//						       FDX			BDX
		// Since a central difference Ix approximation is not possible.
		// This reaches all four corners.
		for (j = 0; j <= iMaxY; j++) {
		  ppfGradient[0][j] = FDX(ppfData,0,j);
		  ppfGradient[iMaxX][j] = BDX(ppfData,iMaxX,j);
        }

        // On all internal pixels we can calculate a central difference approx.
        for (i=1;i < iMaxX; i++) {
		  for(j=1; j < iMaxY; j++) {
		    ppfGradient[i][j] = CDX(ppfData,i,j);
		  }
		}
		// Along edges of image where j == 0 and j == iMaxY (left, right col).
        // we can calculate CDX, but not in the corners.
		for (i = 1; i < iMaxX; i++) {
		  ppfGradient[i][0] = CDX(ppfData,i,0);
          ppfGradient[i][iMaxY] = CDX(ppfData,i,iMaxY);
        }
    }
  return 0; // Success
}


int iGradX(float **ppfData, float **ppfGradient, int numRowX, int numColY)
/* 
 * PURPOSE : Calculate the partial derivative in Y direction.
 * 			X - increases across the rows.
 * 			Y - increases across the columns.
 * INPUTS :
 * ppfData     : The image data.
 * ppfGradient : The gradient will be stored in this array
 * numRowX	   : The number of points in the X direction.
 * numColY	   : The number of points in the Y direction.
 *
 * OUTPUTS :
 * 		Iy is approximated.
 *
 * METHOD : At each point in the image a central difference approx.
 * 	    to the first derivative in the y direction is made.
 * 	    At the boundaries a forward difference or a backward difference
 * 	    approximation is made.
 */
{
    register int i; register int j; register int iMaxX; register int iMaxY;

    iMaxX = numRowX - 1;
    iMaxY = numColY - 1;
    assert(iMaxX >= 0 && iMaxY >= 0);

    if (iMaxY == 0) { /* Iy == 0 everywhere for one column data */
	  for (i = 0; i <= iMaxX; i++) {
		ppfGradient[i][0] = 0;
	  }
    } else { /* iMaxX can still be zero - but this is OK */
		// Along edges of image where j == 0 and j == iMaxY (left and right col)
		//						       FDY			BDY
		// Since a central difference Iy approximation is not possible.
		// This reaches all four corners.
		for (i = 0; i <= iMaxX; i++) {
		  ppfGradient[i][0] = FDY(ppfData,i,0);
		  ppfGradient[i][iMaxY] = BDY(ppfData,i,iMaxY);
        }

        // On all internal pixels we can calculate a central difference approx.
        for (i=1;i < iMaxX; i++) {
		  for(j=1; j < iMaxY; j++) {
		    ppfGradient[i][j] = CDY(ppfData,i,j);
		  }
		}
		// Along edges of image where i == 0 and i == iMaxX (top, bottom rows).
        // we can calculate CDY, but not in the corners.
		for (j = 1; j < iMaxY; j++) {
		  ppfGradient[0][j] = CDY(ppfData,0,j);
          ppfGradient[iMaxX][j] = CDY(ppfData,iMaxX,j);
        }
    }
  return 0; // Success
}

int GradX(float **ppfData, float **ppfGradient, int max1, int max2) {
  int  M1Row, M2Row, M1Col, M2Col,i;
  float **filter;
  M1Row = -1;M2Row = 1;
  M1Col = -1;M2Col = 1;
  filter = (float**)malloc(sizeof(float*)*(M2Row-M1Row+1));
  for (i=0;i<M2Row-M1Row+1;i++) {
    filter[i] = (float*)malloc(sizeof(float)*(M2Row-M1Row+1));
  }
  
  filter[0][0] = -sqrt(2);filter[0][1] = 0;filter[0][2] = sqrt(2);
  filter[1][0] = -1      ;filter[1][1] = 0;filter[1][2] = 1;
  filter[2][0] = -sqrt(2);filter[2][1] = 0;filter[2][2] = sqrt(2);

  Convolution(ppfData, filter, M1Row, M2Row, M1Col, M2Col, ppfGradient,max1,max2);
  free(filter);
}

int GradY(float **ppfData, float **ppfGradient, int max1, int max2) {
  int  M1Row, M2Row, M1Col, M2Col,i;
  float **filter;
  M1Row = -1;M2Row = 1;
  M1Col = -1;M2Col = 1;
  filter = (float**)malloc(sizeof(float*)*(M2Row-M1Row+1));
  for (i=0;i<M2Row-M1Row+1;i++) {
    filter[i] = (float*)malloc(sizeof(float)*(M2Row-M1Row+1));
  }
  
  filter[0][0] = -sqrt(2);filter[0][1] = -1;filter[0][2] = -sqrt(2);
  filter[1][0] = 0      ;filter[1][1] = 0;filter[1][2] = 0;
  filter[2][0] = sqrt(2);filter[2][1] = 1;filter[2][2] = sqrt(2);

  Convolution(ppfData, filter, M1Row, M2Row, M1Col, M2Col, ppfGradient,max1,max2);
  free(filter);
}

int Convolution(float **image, // image to be convolved
	float **filter,  // Array of filter coefficients
	int M1Row, int M2Row, // filter left and right limits
	int M1Col, int M2Col, // filter bottom and top limits
        float **result, int max1, int max2)
// On exit the data in result is the convolution of the filter
// with the image data
// No extensions to the image data are made i.e. some of the overlap
// of the convolution is thrown away
{
  int n; int m; int k; int l;
  float sum = 0.0;

  for (n = 0; n < max1; n++) {
    for (m = 0; m < max2; m++) {
      sum = 0.0;
      for (k = M1Row; k <= M2Row; k++) {
        if (n + k < 0 || n + k > max1 - 1) continue;
        for (l = M1Col; l <= M2Col; l++) {
          if (m + l < 0 || m + l > max2 - 1) continue;
          sum += filter[k - M1Row][l - M1Col] * image[n + k][m + l];
        }
      }
      result[n][m] = sum;
    }
  }
  return 0;
}

int normalize(float** gradientx, float** gradienty, int max1,int max2) {
  int i,j;
  float norma;

  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      norma = sqrt(gradientx[i][j]*gradientx[i][j] + gradienty[i][j]*gradienty[i][j]);
      if (norma != 0) {
	gradientx[i][j] = gradientx[i][j] /norma;
	gradienty[i][j] = gradienty[i][j] /norma;
      }
    }
  }
  return 0;  
}

int compute_corners(unsigned short *input, int max1,int max2) {
  int i,j,k,x,y,start,end,mapindex,newmapindex,flag,num_corners = 0;
  int in,ip,corner[10];
  int count,count0,count1;
  float ax,ay,bx,by;
  unsigned char* aux;
  struct list {
    int num_elem;
    int *elem;
  };
  struct list mylist,mylist2;
  int max_num_list = 2000, num_vec = 10;
  int *ideal;
  float* coseno_ideal;
  float coseno[num_vec+1];

  aux = (unsigned char*)calloc(max1*max2,sizeof(unsigned char));
  mylist.elem = (int *)malloc(max_num_list*sizeof(int));
  mylist.num_elem=0;
  mylist2.elem = (int *)malloc(max_num_list*sizeof(int));
  mylist2.num_elem=0;

  i = 0;
  while (input[i] != 1) {
    i++;
  }
  start = i;
  aux[start] = 1;
  mylist.elem[mylist.num_elem] = start;
  mylist.num_elem++;
  
  mylist2.elem[mylist2.num_elem] = start;
  mylist2.num_elem++;
  while (mylist2.num_elem !=0) {  
    mapindex = mylist2.elem[mylist2.num_elem-1];
    mylist2.num_elem--;
    for (x=-1;x<2;x++) {      
      for (y=-1;y<2;y++) {
	if (x==0 && y ==0) continue;
	newmapindex = mapindex + y*max2 + x;
	if (input[newmapindex] == 1 && aux[newmapindex] == 0) {
	  aux[newmapindex] = 1;	  	 
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  mylist2.elem[mylist2.num_elem] = newmapindex;
	  mylist2.num_elem++;
	  if (mylist.num_elem >  max_num_list) {
	    printf("Error mylist.num_elem > %d \n",max_num_list);
	    return 1;
	  }
	}
      }
    }
  }
  end = mylist.elem[mylist.num_elem-1];

  printf("num_elem en la curva %d\n",mylist.num_elem-1);
  /* printf("start x %d y %d\n",maptox(start,max2),maptoy(start,max2));
  printf("end x %d y %d\n",maptox(end,max2),maptoy(end,max2));
  printf("input[end] %d input[start] %d\n",input[end],input[start]);*/

  /* tengo que crear una lista de puntos pertenecientes a la curva 
     correspondientes a la otra parte*/
  mylist2.elem[mylist2.num_elem] = end;
  mylist2.num_elem++;
  while (mylist2.num_elem !=0) {  
    mapindex = mylist2.elem[mylist2.num_elem-1];
    mylist2.num_elem--;
    for (x=-1;x<2;x++) {
      for (y=-1;y<2;y++) {
	if (x==0 && y ==0) continue;
	newmapindex = mapindex + y*max2 + x;
	if (input[newmapindex] == 0 && aux[newmapindex] == 0) {
	  aux[newmapindex] = 1;    
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  mylist2.elem[mylist2.num_elem] = newmapindex;
	  mylist2.num_elem++;
	  if (mylist.num_elem >  max_num_list) {
	    printf("Error mylist.num_elem > %d \n",max_num_list);
	    return 1;
	  }
	}
      }
    }
  }  

  printf("num_elem en la curva %d\n",mylist.num_elem-1);
  coseno_ideal = (float*)malloc(sizeof(float)*mylist.num_elem);
  ideal = (int*)malloc(sizeof(int)*mylist.num_elem);
  /* tengo que calcular (ax,ay), (bx,by) en una vecindad de num_vec puntos*/ 
  for (i=0;i<mylist.num_elem;i++) {
    for (k = num_vec; k > 0; k--) {
      ip = (i+k) % mylist.num_elem;
      in = (i-k+mylist.num_elem) % mylist.num_elem;
      ax = maptox(mylist.elem[i],max2) - maptox(mylist.elem[ip],max2);
      ay = maptoy(mylist.elem[i],max2) - maptoy(mylist.elem[ip],max2);
      bx = maptox(mylist.elem[i],max2) - maptox(mylist.elem[in],max2);
      by = maptoy(mylist.elem[i],max2) - maptoy(mylist.elem[in],max2);
      coseno[k] = (float)(ax * bx + ay * by) / (sqrt(ax*ax + ay*ay) * sqrt(bx*bx + by*by));      
      if (k < num_vec) {
	if (coseno[k] >= coseno[k+1]) {
	  ideal[i] = k+1;
	  coseno_ideal[i] = (float)coseno[k+1];
	  k = 0;
	}
      } else {
	ideal[i] = k;
	coseno_ideal[i] = coseno[k];
      }
    }
  }
  /* Necesito otra pasada para encontrar los puntos de max curvatura*/
  for (i=0;i<mylist.num_elem;i++) {
    flag = 1;
    for (k=-num_vec/2;k<=num_vec/2;k++) {
      if (k==0) continue;
      in = (i+k+mylist.num_elem) % mylist.num_elem;	
      if (coseno_ideal[in] >= coseno_ideal[i]) {
	flag = 0;
      }
    }
    if (flag == 1 && coseno_ideal[i] >0) {
      printf("found corner at x %d y %d coseno_ideal[i] %f\n",maptox(mylist.elem[i],max2),maptoy(mylist.elem[i],max2),coseno_ideal[i]);
      aux[mylist.elem[i]] = 255;
      corner[num_corners] = i;
      num_corners++; 
      if (num_corners > 10) {
	printf("Error, num_corners > %d\n",num_corners);
	return 1;
      }
      input[mylist.elem[i]] = 7;      
    }  
    input[mylist.elem[i]] = i;
  }

  return 0;
  k = (corner[0] + 1) % mylist.num_elem;
  count0=0;
  count1=0;
  for (i=0;i<mylist.num_elem;i++) {
    mapindex = mylist.elem[k];
    if (aux[mapindex] == 255) {
      printf("count0 %d count1 %d at %d %d\n",count0,count1,maptox(mapindex,max2),maptoy(mapindex,max2));
      count0=0;
      count1=0;
    }
    if (input[mapindex] == 0) {	
      count0++;
    }
    if (input[mapindex] == 1) {	
      count1++;
    }
    k = (k + 1) % mylist.num_elem;
  }
  printf("count0 %d count1 %d\n",count0,count1);

  
  free(coseno_ideal);
  free(ideal);
  free(aux);
  free(mylist.elem);
  free(mylist2.elem);
  
  return 0;
}

int initialize_list(struct fifo *mylist, int max_num_elem) { 
  mylist->last=-1;
  mylist->first=0;
  mylist->index_elem=(int*)malloc(sizeof(int)*max_num_elem);
}

int get_from_list(struct fifo *mylist) {
  int elem;
  if (num_elem_in_list(*mylist) >0) {    
    elem = mylist->index_elem[mylist->first];
    mylist->first++;
    return elem;
  } else {
    printf("Error trying to get element no existent in list\n");
    return -1;
  }
}

int put_in_list(struct fifo *mylist, int element) { 
  mylist->last++;
  mylist->index_elem[mylist->last] = element;
}

int num_elem_in_list(struct fifo mylist) { 
  return mylist.last - mylist.first + 1;
}

int reset_list(struct fifo *mylist) {
  mylist->last=-1;
  mylist->first=0;
}

int delete_list(struct fifo *mylist) { 
  mylist->last=-1;
  mylist->first=0;
  free(mylist->index_elem);
}

int new_compute_corners(unsigned short *input, int max1,int max2) {
  int i,j,k,x,y,start,end,mapindex,newmapindex,flag,num_corners = 0;
  int in,ip,tramo,corner[10];
  int count,count0[10],count1[10];
  float ax,ay,bx,by;
  unsigned char* aux;
  struct list {
    int num_elem;
    int *elem; 
  };
  struct list mylist;
  struct fifo mylist2;
  int max_num_list = 2000, num_vec = 10;
  int *ideal;
  float* coseno_ideal;
  float coseno[num_vec+1];

  aux = (unsigned char*)calloc(max1*max2,sizeof(unsigned char));
  mylist.elem = (int *)malloc(max_num_list*sizeof(int));
  mylist.num_elem=0;
  initialize_list(&mylist2,max_num_list);

  i = 0;
  while (input[i] != 1 && i<max1*max2) {
    i++;
  }
  start = i;
  if (start == max1*max2) {
    /* printf("desapareciendo del mapa \n");*/
    free(mylist.elem);
    free(aux);
    delete_list(&mylist2);
    return 0;
  }
  aux[start] = 1;
  mylist.elem[mylist.num_elem] = start;
  mylist.num_elem++;

  put_in_list(&mylist2,start);
  
  while (num_elem_in_list(mylist2) !=0) {  
    mapindex = get_from_list(&mylist2);
    for (x=-1;x<2;x++) {      
      for (y=-1;y<2;y++) {
	if (x==0 && y ==0) continue;
	newmapindex = mapindex + y*max2 + x;
	if (input[newmapindex] == 1 && aux[newmapindex] == 0) {
	  aux[newmapindex] = 1;	  		  
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  put_in_list(&mylist2,newmapindex);
	  if (mylist.num_elem >  max_num_list) {
	    printf("Error mylist.num_elem > %d \n",max_num_list);
	    return 1;
	  }
	}
      }
    }
  }
  end = mylist.elem[mylist.num_elem-1];

  /* printf("num_elem en la curva %d\n",mylist.num_elem);*/
  /* printf("start x %d y %d\n",maptox(start,max2),maptoy(start,max2));
  printf("end x %d y %d\n",maptox(end,max2),maptoy(end,max2));
  printf("input[end] %d input[start] %d\n",input[end],input[start]);*/

  /* tengo que crear una lista de puntos pertenecientes a la curva 
     correspondientes a la otra parte*/
  reset_list(&mylist2);
  put_in_list(&mylist2,end);
  while (num_elem_in_list(mylist2) !=0) {  
    mapindex = get_from_list(&mylist2);
    for (x=-1;x<2;x++) {
      for (y=-1;y<2;y++) {
	if (x==0 && y ==0) continue;
	newmapindex = mapindex + y*max2 + x;
	if (input[newmapindex] == 0 && aux[newmapindex] == 0) {
	  aux[newmapindex] = 1;
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  put_in_list(&mylist2,newmapindex);
	  if (mylist.num_elem >  max_num_list) {
	    printf("Error mylist.num_elem > %d \n",max_num_list);
	    return 1;
	  }
	}
      }
    }
  }  

  printf("num_elem en la curva %d\n",mylist.num_elem);
  coseno_ideal = (float*)malloc(sizeof(float)*mylist.num_elem);
  ideal = (int*)malloc(sizeof(int)*mylist.num_elem);
  /* tengo que calcular (ax,ay), (bx,by) en una vecindad de num_vec puntos*/ 
  for (i=0;i<mylist.num_elem;i++) {
    for (k = num_vec; k > 0; k--) {
      ip = (i+k) % mylist.num_elem;
      in = (i-k+mylist.num_elem) % mylist.num_elem;
      ax = maptox(mylist.elem[i],max2) - maptox(mylist.elem[ip],max2);
      ay = maptoy(mylist.elem[i],max2) - maptoy(mylist.elem[ip],max2);
      bx = maptox(mylist.elem[i],max2) - maptox(mylist.elem[in],max2);
      by = maptoy(mylist.elem[i],max2) - maptoy(mylist.elem[in],max2);
      coseno[k] = (float)(ax * bx + ay * by) / (sqrt(ax*ax + ay*ay) * sqrt(bx*bx + by*by));      
      if (k < num_vec) {
	if (coseno[k] >= coseno[k+1]) {
	  ideal[i] = k+1;
	  coseno_ideal[i] = (float)coseno[k+1];
	  k = 0;
	}
      } else {
	ideal[i] = k;
	coseno_ideal[i] = coseno[k];
      }
    }
  }
  /* Necesito otra pasada para encontrar los puntos de max curvatura*/
  for (i=0;i<mylist.num_elem;i++) {
    flag = 1;
    for (k=-num_vec/2;k<=num_vec/2;k++) {
      if (k==0) continue;
      in = (i+k+mylist.num_elem) % mylist.num_elem;	
      if (coseno_ideal[in] >= coseno_ideal[i]) {
	flag = 0;
      }
    }
    if (flag == 1 && coseno_ideal[i] >0) {
      printf("found corner at x %d y %d coseno_ideal[i] %f\n",maptox(mylist.elem[i],max2),maptoy(mylist.elem[i],max2),coseno_ideal[i]);
      aux[mylist.elem[i]] = 255;
      corner[num_corners] = i;
      num_corners++; 
      if (num_corners > 10) {
	printf("Error, num_corners > %d\n",num_corners);
	return 1;
      }
      input[mylist.elem[i]] = 7;      
    }  
    /* input[mylist.elem[i]] = i;*/
  }

  if (num_corners < 2) {
     free(coseno_ideal);
     free(ideal);
     free(aux);
     free(mylist.elem);
     delete_list(&mylist2);
     return 0;
  }

  for (i=0;i<10;i++) {
    count0[i] = 0;
    count1[i] = 0;
  }
  k = (corner[0] + 1) % mylist.num_elem;
  tramo = 0;
  for (i=0;i<mylist.num_elem;i++) {
    mapindex = mylist.elem[k];
    if (aux[mapindex] == 255) {
      tramo++;     
    }
    if (input[mapindex] == 0) {	
      count0[tramo]++;
    }
    if (input[mapindex] == 1) {	
      count1[tramo]++;
    }
    k = (k + 1) % mylist.num_elem;
  }
  /* printf("count0[0] %d count1[0] %d\n",count0[0],count1[0]);
     printf("count0[1] %d count1[1] %d\n",count0[1],count1[1]);*/

  if (num_corners == 2) {
    if (count0[0] > count0[1]) {
      count0[0] = 0;
      count0[1] = 1;
    } else {
      count0[0] = 1;
      count0[1] = 0;
    }
    k = (corner[0] + 1) % mylist.num_elem;
    tramo = 0;
    for (i=0;i<mylist.num_elem;i++) {
      mapindex = mylist.elem[k];
      if (aux[mapindex] == 255) {
	tramo++;     
      }      	
      input[mapindex] = count0[tramo];      
      k = (k + 1) % mylist.num_elem;
    }
  } else {
    k = (corner[0] + 1) % mylist.num_elem;
    tramo = 0;
    for (i=0;i<mylist.num_elem;i++) {
      mapindex = mylist.elem[k];
      if (aux[mapindex] == 255) {
	tramo++;     
      }
      if (count0[tramo] > count1[tramo]) {	
	input[mapindex] = 0;
      }
      if (count0[tramo] < count1[tramo]) {	
	input[mapindex] = 1;
      }
      k = (k + 1) % mylist.num_elem;
    }
  }

  free(coseno_ideal);
  free(ideal);
  free(aux);
  free(mylist.elem);
  delete_list(&mylist2);
  
  return 0;
}

int sumar_l1l2(float* input1, float* input2, float* output, int totdim) {
  int i;

  for (i=0;i<totdim;i++) {
    output[i] = input1[i] + input2[i] - 1;
  }

  return 0;

}
