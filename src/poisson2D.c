/*(c) Ruben Cardenes Almeida, Valladolid, 8/11/2007 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "poisson2D.h"

/* Input values
255, o 3: Outside domain
1 Exterior boundary 
0 Interior boundary
2 Inside domain 
*/ 
int poisson2D(unsigned char* input,int max1, int max2, float** output, int iterations, float lambda, int reverse, float h, int cortex_label) {
  int i,j,l;
  int sum = 0;
  /* Initialize domain, inside=0, and boundaries values*/
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (input[sum] != 0) {
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
  
  /* Solve Poisson */
  for (l=0;l<iterations;l++) {
    sum = 0;
    for (i=0;i<max1;i++) {
      for (j=0;j<max2;j++) {
	if (input[sum] != 0
	    && i != 0 && i != max1-1 && j!=0 && j != max2-1) {
	  output[i][j] = 0.25 *(output[i-1][j] + output[i+1][j] + output[i][j-1] + output[i][j+1] - h*h);
	  /*output[i][j] = output[i][j]+(lambda+1)*(0.25 *(output[i-1][j] + output[i+1][j] + output[i][j-1] + output[i][j+1]) - output[i][j]);*/
	}
	sum++;
      }
    }
  }
 
  return 0;
}
  

int maximos_locales2D(float** in, unsigned char* out, int max1, int max2) {
  int i,j,x,y,k=0,sum;
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (i == 0 || j == 0 || i == max1-1 || j == max2-1) continue;
      sum = 0;
      for (x=-1;x<2;x++) {
	for (y=-1;y<2;y++) {
	  if (x ==0 && y ==0) continue;
	  if (in[i+x][j+y] < in[i][j]) {
	    sum++;	
	  }     
	}
      }
      if (sum == 8) {
	out[i*max2+j] = (int)in[i][j];
      } else {
	out[i*max2+j] = 0;
      }
      k++;
    }    
  }

  return 0;
}

int minimos_locales2D(float** in, unsigned char* out, int max1, int max2) {
  int i,j,x,y,k=0,sum;
  for (i=0;i<max1;i++) {
    for (j=0;j<max2;j++) {
      if (i == 0 || j == 0 || i == max1-1 || j == max2-1) continue;
      sum = 0;
      for (x=-1;x<2;x++) {
	for (y=-1;y<2;y++) {
	  if (x ==0 && y ==0) continue;
	  if (in[i+x][j+y] > in[i][j]) {
	    sum++;	
	  }     
	}
      }
      if (sum == 8) {
	out[i*max2+j] = (int)in[i][j];
      } else {
	out[i*max2+j] = 0;
      }
      k++;
    }    
  }

  return 0;
}



