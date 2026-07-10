/*(c) Ruben Cardenes Almeida, Boston, 22/3/2004 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "laplace3D.h"

int maptox3d(int mapindex,int max1,int max2) {
  return (mapindex % (max1*max2)) % max2;
}

int maptoy3d(int mapindex,int max1,int max2) {
  return ((mapindex % (max1*max2)) - ((mapindex % (max1*max2)) % max2) ) / max2;
}

int maptoz3d(int mapindex,int max1,int max2) {
  return (mapindex - (mapindex % (max1*max2))) / (max1*max2);
}

/* Input values
255 Outside domain
1 Exterior boundary 
0 Interior boundary
2 Inside domain 
*/ 
int laplace3D(unsigned char* input,int max1, int max2, int max3, float*** output, int iterations, float lambda) {
  int i,j,k,l;
  int sum = 0;
  /* Initialize domain, inside=0, and boundaries values*/
  for (k=0;k<max3;k++) {
    for (i=0;i<max1;i++) {
      for (j=0;j<max2;j++) {
	if (input[sum] == 2) {
	  output[k][i][j] = 0;
	} else if (input[sum] == 3) {
	  output[k][i][j] = -1;
	} else {
	  output[k][i][j] = input[sum];
	}
	sum++;
      }
    }
  }
  
  /* Solve Laplacian */
  for (l=0;l<iterations;l++) {
    sum = 0;
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	for (j=0;j<max2;j++) {
	  if (input[sum] == 2
	      && i != 0 && i != max1-1 && j!=0 && j != max2-1 && k != 0 && k != max3-1) {
	    /*output[k][i][j] = (output[k-1][i][j] + output[k+1][i][j] + output[k][i-1][j] + output[k][i+1][j] + output[k][i][j-1] + output[k][i][j+1])/6;*/
	    output[k][i][j] = output[k][i][j]+(lambda+1)*(0.16666667 *(output[k-1][i][j] + output[k+1][i][j] + output[k][i-1][j] + output[k][i+1][j] + output[k][i][j-1] + output[k][i][j+1]) - output[k][i][j]);
	  }
	  sum++;
	}
      }
    }
  }
 
  return 0;
}

int laplace3D_voxelsize(unsigned char* input,int max1, int max2, int max3, float*** output, int iterations,  float hx, float hy, float hz, float lambda) {
  int i,j,k,l;
  int sum = 0;
  /* Initialize domain, inside=0, and boundaries values*/
  for (k=0;k<max3;k++) {
    for (i=0;i<max1;i++) {
      for (j=0;j<max2;j++) {
	if (input[sum] == 2) {
	  output[k][i][j] = 0;
	} else if (input[sum] == 3) {
	  output[k][i][j] = -1;
	} else {
	  output[k][i][j] = input[sum];
	}
	sum++;
      }
    }
  }
  
  /* Solve Laplacian */
  for (l=0;l<iterations;l++) {
    sum = 0;
    for (k=0;k<max3;k++) {
      for (i=0;i<max1;i++) {
	for (j=0;j<max2;j++) {
	  if (input[sum] == 2
	      && i != 0 && i != max1-1 && j!=0 && j != max2-1 && k != 0 && k != max3-1) {
	    output[k][i][j] = 0.5 *( (output[k-1][i][j] + output[k+1][i][j])/(hz*hz) + (output[k][i-1][j] + output[k][i+1][j])/(hy*hy) + (output[k][i][j-1] + output[k][i][j+1])/(hx*hx))*(hx*hx*hy*hy*hz*hz)/(hx*hx*hy*hy + hy*hy*hz*hz + hx*hx*hz*hz);
	    /* output[k][i][j] = output[k][i][j]+(lambda+1)*(0.5 *( (output[k-1][i][j] + output[k+1][i][j])/(hz*hz) + (output[k][i-1][j] + output[k][i+1][j])/(hy*hy) + (output[k][i][j-1] + output[k][i][j+1])/(hx*hx))*(hx*hx*hy*hy*hz*hz)/(hx*hx*hy*hy + hy*hy*hz*hz + hx*hx*hz*hz) - output[k][i][j]);*/
	  }
	  sum++;
	}
      }
    }
  }
 
  return 0;
}

int EdgeDetect3D_knee(unsigned char *domain, int max1, int max2, int max3) {
  int x,y,z,i;
  i = 0;
 
  for (z=0;z<max3;z++) {
    for(y=0; y<max2; y++) {
      for(x=0; x<max1; x++) {
	if ((x==0)||(y==0)||(z==0)||(x==max1-1)||(y==max2-1)||(z==max3-1)) {	 
	  /* domain[i]=255;*/	 
	} else
	  if ( (domain[i]==2)&&
	       ((domain[i+1]==255)||
		(domain[i-1]==255)||

		(domain[i+max1]==255)||
		(domain[i-max1]==255)||

		(domain[i+max1*max2]==255)||
		(domain[i-max1*max2]==255))) {
	    
	    domain[i]=1;	    	    
	  }
	/*else {	   	   
	  domain[i]=255;
	  }*/
	i++;
      }
    }
  }
  return 0;
}

int EdgeDetect3D(unsigned char *domain, int max1, int max2, int max3) {
  int x,y,z,i;
  i = 0;
 
  for (z=0;z<max3;z++) {
    for(y=0; y<max2; y++) {
      for(x=0; x<max1; x++) {
	if ((x==0)||(y==0)||(z==0)||(x==max1-1)||(y==max2-1)||(z==max3-1)) {	 
	  /* domain[i]=255;*/	 
	} else
	  if ( (domain[i]!=0)&&
	       ((domain[i+1]==0)||
		(domain[i-1]==0)||

		(domain[i+max1]==0)||
		(domain[i-max1]==0)||

		(domain[i+max1*max2]==0)||
		(domain[i-max1*max2]==0))) {
	  
	    domain[i]=1;	    	    
	  }
	/*else {	   	   
	  domain[i]=255;
	  }*/
	i++;
      }
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

int relabel_ushort(unsigned short* data, int totdim, unsigned char oldlabel, unsigned char newlabel) {
  int i;
  for (i=0;i<totdim;i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
} 

int relabel_uchar(unsigned char* data, int totdim, unsigned char oldlabel, unsigned char newlabel) {
  int i;
  for (i=0;i<totdim;i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
}

int relabel_float(float* data, int totdim, float oldlabel, float newlabel) {
  int i;
  for (i=0;i<totdim;i++) {
    if (data[i] == oldlabel) {
      data[i] = newlabel;
    }
  }
  return 0;
}

int floodfill3D(unsigned char *domain, int startindex, unsigned short oldlabel, unsigned char newlabel,int max1, int max2,int max3) {
  int x,y,z,newmapindex,mapindex,xnew,ynew,znew;
  struct list {
    int num_elem;
    int *elem;
  } mylist;
  int max_num_elem_mylist = max1*max2*max3;

  mylist.elem = (int*)malloc(sizeof(int)*max_num_elem_mylist);
  mylist.elem[0] = startindex;
  mylist.num_elem = 1;  

  while (mylist.num_elem != 0) {
    mapindex = mylist.elem[mylist.num_elem-1];
    mylist.num_elem--;
    for (z=-1;z<2;z=z+2) {
      x=0;y=0;     
      newmapindex = mapindex + max1*max2*z;
      xnew = maptox3d(newmapindex,max1,max2);
      ynew = maptoy3d(newmapindex,max1,max2);
      znew = maptoz3d(newmapindex,max1,max2);
      
      if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 && znew >=0 && znew < max3) {
	if (domain[newmapindex] == oldlabel) {
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  if (mylist.num_elem > max_num_elem_mylist) {
	    printf("Error, exceeded num elem in mylist: %d, in floodfill3D \n", max_num_elem_mylist);
	    return 1;
	  }
	  domain[newmapindex] = newlabel;
	}
      }
    }
    for (x=-1;x<2;x++) {
      z=0;y=0;     
      newmapindex = mapindex + x;
      xnew = maptox3d(newmapindex,max1,max2);
      ynew = maptoy3d(newmapindex,max1,max2);
      znew = maptoz3d(newmapindex,max1,max2);
      
      if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 && znew >=0 && znew < max3) {
	if (domain[newmapindex] == oldlabel) {
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  if (mylist.num_elem > max_num_elem_mylist) {
	    printf("Error, exceeded num elem in mylist: %d, in floodfill3D \n", max_num_elem_mylist);
	    return 1;
	  }
	  domain[newmapindex] = newlabel;
	}
      }
    }
    for (y=-1;y<2;y=y+2) {  
      x=0;z=0;     
      newmapindex = mapindex + max2*y;
      xnew = maptox3d(newmapindex,max1,max2);
      ynew = maptoy3d(newmapindex,max1,max2);
      znew = maptoz3d(newmapindex,max1,max2);
      
      if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 && znew >=0 && znew < max3) {
	if (domain[newmapindex] == oldlabel) {
	  mylist.elem[mylist.num_elem] = newmapindex;
	  mylist.num_elem++;
	  if (mylist.num_elem > max_num_elem_mylist) {
	    printf("Error, exceeded num elem in mylist: %d, in floodfill3D \n", max_num_elem_mylist);
	    return 1;
	  }
	  domain[newmapindex] = newlabel;
	}
      }
    }
  }

  free(mylist.elem);
  return 0;
}

int RelabelBoundary3D(unsigned char *domain,int max1,int max2,int max3){
  int x,y,z,xr,yr,zr,i,xnew,ynew,znew;
  int newmapindex,mapindex,start,newstart;
  struct list {
    int num_elem;
    int *elem;
  } mylist;
  int max_num_elem_mylist = max1*max2*max3;

  i = 0;
  
  for(z=0; z<max3; z++) {
    for(y=0; y<max2; y++) {
      for(x=0; x<max1; x++) {
	if (domain[i] == 1) {
	  start = i;
	  y=max1;x=max2;z=max3;
	}
	i++;    
      }
    }
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
    for (z=-1;z<2;z++) {
      for (x=-1;x<2;x++) {
	for (y=-1;y<2;y++) {
	  if (x==0 && y ==0 && z==0) continue;
	  newmapindex = mapindex + z*max1*max2 + y*max2 + x;
	  xnew = maptox3d(newmapindex,max1,max2);
	  ynew = maptoy3d(newmapindex,max1,max2);
	  znew = maptoz3d(newmapindex,max1,max2);
	  if (xnew >= 0 && xnew < max2 && ynew >= 0 && ynew < max1 && znew >= 0 && znew < max3 &&
	      domain[newmapindex] == 1) {
	    /* Put new element in mylist*/
	    mylist.elem[mylist.num_elem] = newmapindex;
	    mylist.num_elem++;
	    if (mylist.num_elem > max_num_elem_mylist) {
	      printf("Error, exceeded number of elements in mylist.num_elem ( > %d ), in RelabelBoundary3D\n",max_num_elem_mylist);
	      return 1;
	    }
	    domain[newmapindex] = 0;
	  }
	}
      }
    }
  }
  printf("Doing floodfill3D, domain[0] %d \n",domain[0]);
  /* floodill the exterior region, relabeling to -1 */ 
  floodfill3D(domain,0,domain[0],3,max1,max2,max3);
  free(mylist.elem);
  return 0;
}


int iGradX3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY,int numSlice,float hx)
/* 
 * PURPOSE : Calculate the partial derivative in X direction.
 * 			X - increases across the rows.
 * 			Y - increases across the columns.
 * 			Z - increases across the slices.
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
    register int i; register int j;  register int k;register int iMaxX; register int iMaxY;register int iMaxZ;

    iMaxX = numRowX - 1;
    iMaxY = numColY - 1;
    iMaxZ = numSlice - 1;
    assert(iMaxX >= 0 && iMaxY >= 0 && iMaxZ >= 0);

    if (iMaxX == 0) { /* Ix == 0 everywhere for one row data */
      for (k = 0; k <= iMaxZ; k++) {
	for (i = 0; i <= iMaxY; i++) {
	  ppfGradient[k][i][0] = 0;
	}
      }
    } else { /* iMaxY can still be zero - but this is OK */
		// Along edges of image where i == 0 and i == iMaxX (top and bottom row)
		//						       FDX			BDX
		// Since a central difference Ix approximation is not possible.
		// This reaches all four corners.
      for (k = 0; k <= iMaxZ; k++) {
	for (i = 0; i <= iMaxY; i++) {
	  ppfGradient[k][i][0] = FDX3(ppfData,k,i,0,hx);
	  ppfGradient[k][i][iMaxX] = BDX3(ppfData,k,i,iMaxX,hx);
	}
      }

      // On all internal pixels we can calculate a central difference approx.
      for (k=0;k <= iMaxZ; k++) {
	for (i=0;i<=iMaxY;i++) {
	  for(j=1;j<iMaxX;j++) {
	    ppfGradient[k][i][j] = CDX3(ppfData,k,i,j,hx);	   
	  }
	}      
      }
    }
  return 0; // Success
}


int iGradY3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY,int numSlice,float hy)
/* 
 * PURPOSE : Calculate the partial derivative in Y direction.
 * 			X - increases across the columns.
 * 			Y - increases across the rows. 
 * 			Z - increases across the slices.
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
    register int i; register int j;  register int k;register int iMaxX; register int iMaxY;register int iMaxZ;

    iMaxX = numRowX - 1;
    iMaxY = numColY - 1;
    iMaxZ = numSlice - 1;
    assert(iMaxX >= 0 && iMaxY >= 0 && iMaxZ >= 0);

    if (iMaxY == 0) { /* Iy == 0 everywhere for one column data */
      for (k = 0; k <= iMaxZ; k++) {
	for (j = 0; j <= iMaxX; j++) {
	  ppfGradient[k][0][j] = 0;
	}
      }
    } else { /* iMaxX can still be zero - but this is OK */
		// Along edges of image where j == 0 and j == iMaxY (left and right col)
		//						       FDY			BDY
		// Since a central difference Iy approximation is not possible.
		// This reaches all four corners.
      for (k = 0; k <= iMaxZ; k++) {
	for (j = 0; j <= iMaxX; j++) {
	  ppfGradient[k][0][j] = FDY3(ppfData,k,0,j,hy);
	  ppfGradient[k][iMaxY][j] = BDY3(ppfData,k,iMaxY,j,hy);
	}
      }

      // On all internal pixels we can calculate a central difference approx.
      for (k=0;k<=iMaxZ;k++) {
	for (i=1;i<iMaxY;i++) {
	  for(j=0;j<=iMaxX;j++) {
	    ppfGradient[k][i][j] = CDY3(ppfData,k,i,j,hy);
	  }
	}
      }      
    }
  return 0; // Success
}

int iGradZ3D(float ***ppfData, float ***ppfGradient, int numRowX, int numColY,int numSlice,float hz)
/* 
 * PURPOSE : Calculate the partial derivative in Z direction.
 * 			X - increases across the columns.
 * 			Y - increases across the rows. 
 * 			Z - increases across the slices.
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
    register int i; register int j;  register int k;register int iMaxX; register int iMaxY;register int iMaxZ;

    iMaxX = numRowX - 1;
    iMaxY = numColY - 1;
    iMaxZ = numSlice - 1;
    assert(iMaxX >= 0 && iMaxY >= 0 && iMaxZ >= 0);

    if (iMaxZ == 0) { /* Iy == 0 everywhere for one column data */
      for (i = 0; i <= iMaxY; i++) {
	for (j = 0; j <= iMaxX; j++) {
	  ppfGradient[0][i][j] = 0;
	}
      }
    } else { /* iMaxX can still be zero - but this is OK */
		// Along edges of image where j == 0 and j == iMaxY (left and right col)
		//						       FDY			BDY
		// Since a central difference Iy approximation is not possible.
		// This reaches all four corners.
      for (i = 0; i <= iMaxY; i++) {
	for (j = 0; j <= iMaxX; j++) {
	  ppfGradient[0][i][j] = FDZ3(ppfData,0,i,j,hz);
	  ppfGradient[iMaxZ][i][j] = BDZ3(ppfData,iMaxZ,i,j,hz);
	}
      }

      // On all internal pixels we can calculate a central difference approx.
      for (k=1;k<iMaxZ;k++) {
	for (i=0;i<=iMaxY;i++) {
	  for(j=0;j<=iMaxX;j++) {
	    ppfGradient[k][i][j] = CDZ3(ppfData,k,i,j,hz);
	  }
	}
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

int normalize3D(float*** gradientx, float*** gradienty,float*** gradientz,int max1,int max2, int max3) {
  int i,j,k;
  float norma;

  for (k=0;k<max3;k++) {
    for (i=0;i<max1;i++) {
      for (j=0;j<max2;j++) {
	norma = sqrt(gradientx[k][i][j]*gradientx[k][i][j] + gradienty[k][i][j]*gradienty[k][i][j]
		     + gradientz[k][i][j]*gradientz[k][i][j]);
	if (norma != 0) {
	  gradientx[k][i][j] = gradientx[k][i][j] /norma;
	  gradienty[k][i][j] = gradienty[k][i][j] /norma;
	  gradientz[k][i][j] = gradientz[k][i][j] /norma;
	}
      }
    }
  }
  return 0;
  
}

int sumar_l1l2(float* input1, float* input2, float* output, int totdim) {
  int i;

  for (i=0;i<totdim;i++) {
    output[i] = input1[i] + input2[i] - 1;
  }

  return 0;

}
