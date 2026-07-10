/* Copyright (c) Ruben Cardenes Almeida 08/04/2004 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "thickness2D.h"
#include "laplace2D.h"

#define PI 3.1415927

static int numelembucket[NUM_BUCKETS];
extern int numrechazos;
extern int numasignaciones;
extern int asignacionesraras;
extern int numPrototypes;
int highestIndexClass;
int actualDimension; 
int numPrototypesInClass[MAXCLASSNUMBER];
char buffer[2048];
int pdim;

int mapIndex3D(int r,int c,int z, int nr,int nc,int nz)
{
  if (c >= nc) return -1;
  if (c < 0) return -1;
  if (r >= nr) return -1;
  if (r < 0) return -1;
  if (z >= nz) return -1;
  if (z < 0) return -1;
  return c + r * nc + z * nr * nc;
}

int mapIndex2D(int r,int c, int nr,int nc)
{
  if (c >= nc) return -1;
  if (c < 0) return -1;
  if (r >= nr) return -1;
  if (r < 0) return -1;
  return c + r * nc;
}

void print_timing(FILE *fp, struct timeval start, struct timeval end) 
{
  double tuend = 1e-06*(double)end.tv_usec; \
  double tustart = 1e-06*(double)start.tv_usec; \
  double tend = end.tv_sec + tuend;\
  double tstart = start.tv_sec + tustart;\
  fprintf(fp,"Elapsed time: %g\n", (tend - tstart) ); \
}

int maximo(int a,int b) {
  if (a > b) {
    return a;
  } else {
    return b;
  }
}

float distanceYezzi_reverse_orig(float **gradientx,float **gradienty,int newmapindex, int x, int y,float *maps,int max2, float r,float hx, float hy) {
  float distf;

  if (gradientx[y][x] > 0) { 
    if (maps[newmapindex+1] == -1 && fabs(gradientx[y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex+1] == -1) { 
	distf = 1.0;
      } else {
	distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex+1]/hx;	      
      }
    }
  } else {
    if (maps[newmapindex-1] == -1 && fabs(gradientx[y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex-1] == -1) { 
	distf = 1.0;
      } else {
	distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex-1]/hx;
      }
    }
  }
  if (gradienty[y][x] > 0) {
    if (maps[newmapindex+max2] == -1 && fabs(gradienty[y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex+max2] > -1) { 
	distf += fabs(gradienty[y][x])*maps[newmapindex+max2]/hy;
      }
    }
  } else {
    if (maps[newmapindex-max2] == -1 && fabs(gradienty[y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex-max2] > -1) { 
	distf += fabs(gradienty[y][x])*maps[newmapindex-max2]/hy;
      }
    }
  }
  distf = distf/(fabs(gradientx[y][x])/hx + fabs(gradienty[y][x])/hy);
  return distf;
}

float distanceYezzi_orig(float **gradientx,float **gradienty,int newmapindex,int x, int y,float *maps,int max2,float r,float hx, float hy) {
  float distf;
  if (gradientx[y][x] < 0) { 	      
    if (maps[newmapindex+1] == -1 && fabs(gradientx[y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex+1] == -1) {
	distf = 1.0;
      } else {
	distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex+1]/hx;
      }
    }	      
  } else {
    if (maps[newmapindex-1] == -1 && fabs(gradientx[y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex-1] == -1) {
	distf = 1.0;
      } else { 
	distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex-1]/hx;
      }
    }
  }
  if (gradienty[y][x] < 0) {
    if (maps[newmapindex+max2] == -1 && fabs(gradienty[y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex+max2] > -1) {
	distf += fabs(gradienty[y][x])*maps[newmapindex+max2]/hy;
      }
    }
  } else {
    if (maps[newmapindex-max2] == -1 && fabs(gradienty[y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex-max2] > -1) {
	distf += fabs(gradienty[y][x])*maps[newmapindex-max2]/hy;
      }
    }
  }
  distf = distf/(fabs(gradientx[y][x])/hx + fabs(gradienty[y][x])/hy);
  return distf;
}

float distanceYezzi_reverse(float **gradientx,float **gradienty,int newmapindex, int x, int y,float *maps,int max2, float r,float hx, float hy) {
  float distf;

  if (gradientx[y][x] > 0) { 
    if (maps[newmapindex+1] == -1) { 
	  distf = 1.0;
	} else {
	  distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex+1]/hx;	      
	}
  } else {
	if (maps[newmapindex-1] == -1) { 
	  distf = 1.0;
	} else {
	  distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex-1]/hx;
	}
  }
  if (gradienty[y][x] > 0) {
	if (maps[newmapindex+max2] > -1) { 
	  distf += fabs(gradienty[y][x])*maps[newmapindex+max2]/hy;
	}
  } else {
	if (maps[newmapindex-max2] > -1) { 
	  distf += fabs(gradienty[y][x])*maps[newmapindex-max2]/hy;
	}
  }
  distf = distf/(fabs(gradientx[y][x])/hx + fabs(gradienty[y][x])/hy);
  return distf;
}

float distanceYezzi(float **gradientx,float **gradienty,int newmapindex,int x, int y,float *maps,int max2,float r,float hx, float hy) {
  float distf;
  if (gradientx[y][x] < 0) {
	if (maps[newmapindex+1] == -1) {
	  distf = 1.0;
    } else {
	  distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex+1]/hx;
	}
  } else {
	if (maps[newmapindex-1] == -1) {
	  distf = 1.0;
	} else { 
	  distf = 1.0 + fabs(gradientx[y][x])*maps[newmapindex-1]/hx;
	}
  }
  if (gradienty[y][x] < 0) {
	if (maps[newmapindex+max2] > -1) {
	  distf += fabs(gradienty[y][x])*maps[newmapindex+max2]/hy;
	}
  } else {
	if (maps[newmapindex-max2] > -1) {
	  distf += fabs(gradienty[y][x])*maps[newmapindex-max2]/hy;
	}
  }
  distf = distf/(fabs(gradientx[y][x])/hx + fabs(gradienty[y][x])/hy);
  return distf;
}


int thickness2DYezzi(unsigned char* prototypes,int max1, int max2, float *maps, float** laplacefield,float** gradientx, float** gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug) {
  int i,j,x,y,xr,yr,mapindex,newmapindex,aux,d,l,flag;
  int *index_aux;
  float distf,r;
  struct list {
    int num_elem;
    int *elem;
  } list1; 
  struct list list2,list_aux;
  int max_number_in_list = 50000;
  unsigned char* prot_copia;
  
  prot_copia = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2);

  for (j=0;j<max1*max2;j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list1.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list1.num_elem = 0;

  list2.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list2.num_elem = 0;
  
  list_aux.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list_aux.num_elem = 0;
  
  for (l=0;l<num_it;l++) {
    d = 0;
    if (l == 0) {
      r=0.3;
    } else {
      r=0.04;
    }
	if (debug == 1) {
      printf("iteration %d\n",l);
	}
    for (i=0;i<max1*max2;i++) {
      prot_copia[i] = prototypes[i]; /* to recover the domain after the first iteration */
      if (prot_copia[i] == 0) {
	     if (list1.num_elem >= max_number_in_list) {
	       printf("Error list1.num_elem >= %d\n",max_number_in_list);
	       return 1;
	     }
	     list1.elem[list1.num_elem] = i;
	     list1.num_elem++;
	     maps[i]=0;
      }
    }

    while (list1.num_elem != 0 || list2.num_elem != 0) {
      /* printf("num elem list1: %d\n",list1.num_elem);*/
      while (list1.num_elem != 0) {
	    /* Get element from list1 */
	    mapindex = list1.elem[list1.num_elem-1];
	    list1.num_elem--;
	    flag = 0;
	    for (x=-1;x<2;x++) {
	      for (y=-1;y<2;y++) {
	        if (x==0 && y ==0) continue;
	        newmapindex = mapindex + x + y*max2;
	        xr = maptox(newmapindex,max2);
	        yr = maptoy(newmapindex,max2);
	        if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1) continue;
	        if (prot_copia[newmapindex] == label_cortex) {
	          /* Compute new distance */
	          distf = distanceYezzi(gradientx,gradienty,newmapindex,xr,yr,maps,max2,0.00,hx,hy);

	          if (distf > 0) {
		        maps[newmapindex] = distf;
	    	    flag = 1;
    		    /* Put new element in list2*/
	    	    if (list2.num_elem >= max_number_in_list) {
	    	      printf("Error list2.num_elem >= %d\n",max_number_in_list);
	    	      return 1;
	    	    }
	    	    list2.elem[list2.num_elem] = newmapindex;
		        list2.num_elem++;
    		    prot_copia[newmapindex] = 0;
	          }
	        }
    	  }
	    }
	    if (flag == 0) {
	      /* Ponemos en un lista auxiliar el punto que no pudo propagarse*/
	      if (list_aux.num_elem >= max_number_in_list) {
	        printf("Error list_aux.num_elem >= %d\n",max_number_in_list);
	        return 1;
	      }
	      list_aux.elem[list_aux.num_elem] = mapindex;
	      list_aux.num_elem++;
	    }
      }

      for (i=0;i<list_aux.num_elem;i++) {
	  /* while (list_aux.num_elem != 0) {    */
	  /* Get element from list1 */
	  mapindex = list_aux.elem[i];
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      if (x==0 && y ==0) continue;
	      newmapindex = mapindex + x + y*max2;
	      xr = maptox(newmapindex,max2);
	      yr = maptoy(newmapindex,max2);
	      if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1) continue;
	      if (prot_copia[newmapindex] == label_cortex) {
	        /* Compute new distance */
	        distf = distanceYezzi(gradientx,gradienty,newmapindex,xr,yr,maps,max2,r,hx,hy);

	        if (distf > 0) {
		  maps[newmapindex] = distf;
		  /* Put new element in list2*/
		  if (list2.num_elem >= max_number_in_list) {
		    printf("Error list2.num_elem >= %d\n",max_number_in_list);
		    return 1;
		  }
		  list2.elem[list2.num_elem] = newmapindex;
		  list2.num_elem++;
		  prot_copia[newmapindex] = 0;
	      }
	    }
	  }
	}
      }
      list_aux.num_elem = 0;
      /* printf("num elem list2: %d\n",list2.num_elem);*/
      d++;
      /* swap(list1.elem,list2.elem); */
      
      index_aux = list1.elem;
      list1.elem = list2.elem;
      list2.elem = index_aux;
      aux = list1.num_elem;
      list1.num_elem = list2.num_elem;
      list2.num_elem = aux;
    }
  }
  printf("dmax = %d\n",d);

  free(list1.elem);
  free(list2.elem);
  free(list_aux.elem);
  free(prot_copia);
  return 0; /* success */
}

int thickness2DYezzi_reverse(unsigned char* prototypes,int max1, int max2, float *maps, float** laplacefield,float** gradientx, float** gradienty, int num_it, float hx, float hy, unsigned char label_cortex, int debug) {
  int i,j,x,y,xr,yr,mapindex,newmapindex,aux,d,flag,l;
  int *index_aux;
  float distf,r;
  struct list {
    int num_elem;
    int *elem;
  } list1; 
  struct list list2,list_aux;
  int max_number_in_list = 50000;
  unsigned char* prot_copia;
  
  prot_copia = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2);

  for (j=0;j<max1*max2;j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list1.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list1.num_elem = 0;

  list2.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list2.num_elem = 0;
  
  list_aux.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list_aux.num_elem = 0;
 
  for (l=0;l<num_it;l++) {
    d = 0;
    if (l == 0) {
      r=0.08;
    } else {
      r=0.00;
    }
	if (debug == 1) {
      printf("iteration %d\n",l);
	}
    for (i=0;i<max1*max2;i++) {
      prot_copia[i] = prototypes[i]; /* to recover the domain after the first iteration */
      if (prot_copia[i] == 1) {
	if (list1.num_elem >= max_number_in_list) {
	  printf("Error list1.num_elem >= %d\n",max_number_in_list);
	  return 1;
	}
	list1.elem[list1.num_elem] = i;
	list1.num_elem++;
	maps[i]=0;
      }
    }

    d = 0;
    while (list1.num_elem != 0 || list2.num_elem != 0) {
      /* printf("num elem list1: %d\n",list1.num_elem);*/
      while (list1.num_elem != 0) {
	/* Get element from list1 */
	mapindex = list1.elem[list1.num_elem-1];
	list1.num_elem--;
	flag = 0;
	for (x=-1;x<2;x++) {
	  for (y=-1;y<2;y++) {
	    if (x==0 && y ==0) continue;
	    newmapindex = mapindex + x + y*max2;
	    xr = maptox(newmapindex,max2);
	    yr = maptoy(newmapindex,max2);
	    if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1) continue;
	    if (prot_copia[newmapindex] == label_cortex) {
	      /* Compute new distance */
	      distf = distanceYezzi_reverse(gradientx,gradienty,newmapindex,xr,yr,maps,max2,0.000,hx,hy);

	      if (distf > 0) {
		flag = 1;
		maps[newmapindex] = distf;
		/* Put new element in list2*/
		if (list2.num_elem >= max_number_in_list) {
		  printf("Error list2.num_elem >= %d\n",max_number_in_list);
		  return 1;
		}
		list2.elem[list2.num_elem] = newmapindex;
		list2.num_elem++;
		prot_copia[newmapindex] = 1;
	      }
	    }
	  }
	}
	if (flag == 0) {
	  /* Ponemos en un lista auxiliar el punto que no pudo propagarse*/
	  if (list_aux.num_elem >= max_number_in_list) {
	    printf("Error list_aux.num_elem >= %d\n",max_number_in_list);
	    return 1;
	  }
	  list_aux.elem[list_aux.num_elem] = mapindex;
	  list_aux.num_elem++;
	}
      }

      for (i=0;i<list_aux.num_elem;i++) {
	/* while (list_aux.num_elem != 0) { */
	/* Get element from list1 */
	mapindex = list_aux.elem[i];
	for (x=-1;x<2;x++) {
	  for (y=-1;y<2;y++) {
	    if (x==0 && y ==0) continue;
	    newmapindex = mapindex + x + y*max2;
	    xr = maptox(newmapindex,max2);
	    yr = maptoy(newmapindex,max2);
	    if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1) continue;
	    if (prot_copia[newmapindex] == label_cortex) {
	      /* Compute new distance */
	      distf = distanceYezzi_reverse(gradientx,gradienty,newmapindex,xr,yr,maps,max2,r,hx,hy);

	      if (distf > 0) {
		maps[newmapindex] = distf;
		/* Put new element in list2*/
		if (list2.num_elem >= max_number_in_list) {
		  printf("Error list2.num_elem >= %d\n",max_number_in_list);
		  return 1;
		}
		list2.elem[list2.num_elem] = newmapindex;
		list2.num_elem++;
		prot_copia[newmapindex] = 1;
	      }
	    }
	  }
	}
      }
   
      /* printf("num elem list2: %d, list_aux.num_elem %d\n",list2.num_elem,list_aux.num_elem);*/
      list_aux.num_elem = 0;
      d++;
      /* swap(list1.elem,list2.elem); */
      
      index_aux = list1.elem;
      list1.elem = list2.elem;
      list2.elem = index_aux;
      aux = list1.num_elem;
      list1.num_elem = list2.num_elem;
      list2.num_elem = aux;
    }
  }
  printf("dmax = %d\n",d);

  free(list1.elem);
  free(list2.elem);
  free(list_aux.elem);
  free(prot_copia);
  return 0; /* success */
}

float getnewcoordinates(int *newx,int *newy,int x,int y,int max2, float gradientx,float gradienty,float *angle_error) {
  double theta,newtheta,dy,dx;
  float new_angle_error = 0;
  theta = acos((double)gradientx);
  if (gradienty < 0) {
    theta = 2 * PI - theta;
  }

  newtheta = theta + angle_error[y*max2+x];
  dy = sin(newtheta);
  dx = cos(newtheta);

  if (dy >=0) { /* theta -> (0,PI)*/
    if (fabs(dx) <= 0.3826) { /*theta -> (3PI/8,5PI/8) */
      (*newx) = x; 
      (*newy) = y-1; 
      new_angle_error = newtheta - PI / 2;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - PI / 2;*/
    } 
    if (dx >= 0.9238) { /*theta -> (0,PI/8) */
      (*newx) = x+1;
      (*newy) = y;
      new_angle_error = newtheta;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta;*/
    }
    if (dx <= -0.9238) { /*theta -> (7PI/8,PI) */
      (*newx) = x-1;
      (*newy) = y;
      new_angle_error = newtheta - PI;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - PI;*/
    }
    if (dx <= 0.9238 && dx >= 0.3826) { /*theta -> (PI/8,3PI/8) */
      (*newx) = x+1;
      (*newy) = y-1;
      new_angle_error = newtheta - PI / 4;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - PI / 4; */
    }
    if (dx >= -0.9238 && dx <= -0.3826) { /*theta -> (5PI/8,7PI/8) */
      (*newx) = x-1;
      (*newy) = y-1;
      new_angle_error = newtheta - 3*PI / 4;
     /* angle_error[(*newx)+(*newy)*max2] = newtheta - 3*PI / 4; */
    }
  }
  if (dy <0) {
    if (fabs(dx) <= 0.3826) {
      (*newx) = x;
      (*newy) = y+1;
      new_angle_error = newtheta - 3*PI / 2;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - 3*PI / 2; */
    }
    if (dx >= 0.9238) {
      (*newx) = x+1;
      (*newy) = y;
      new_angle_error = newtheta;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta; */
    }
    if (dx <= -0.9238) {
      (*newx) = x-1;
      (*newy) = y;
      new_angle_error = newtheta - PI;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - PI;*/
    }
    if (dx <= 0.9238 && dx >= 0.3826) {
      (*newx) = x+1;
      (*newy) = y+1;
      new_angle_error = newtheta -  7*PI / 4;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - 7*PI / 4; */
    }
    if (dx >= -0.9238 && dx <= -0.3826) {
      (*newx) = x-1;
      (*newy) = y+1; 
      new_angle_error = newtheta - 5*PI / 4;
      /* angle_error[(*newx)+(*newy)*max2] = newtheta - 5*PI / 4; */
    }
  }

  return new_angle_error;
}

int thickness2Dgradient(unsigned char* prototypes,int max1, int max2, float *input_maps,
			float* maps, float** laplacefield, float** gradientx, float** gradienty) {
  int i,j,x,y,xr,yr,newx,newy,provx,provy,mapindex,newmapindex,provmapindex,aux,d;
  int *index_aux;
  int counter = 0;
  float distf,xf,yf,new_angle_error;
  float* angle_error;
  int numelemaislados = 0;
  struct list {
    int num_elem;
    int *elem;
  } list1; 
  struct list list2,list3;
  int max_number_in_list = 50000;

  angle_error = (float*)malloc(sizeof(float)*max2*max2);
  for (j=0;j<max1*max2;j++) {
    maps[j] = -1;
    angle_error[j] = 0;
  }

  list1.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list1.num_elem = 0;

  list2.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list2.num_elem = 0;

  list3.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list3.num_elem = 0;

  for (i=0;i<max1*max2;i=i+10) {
    if (prototypes[i] == 2 &&
	(prototypes[i+1] == 0 || prototypes[i-1] == 0 ||
	 prototypes[i+max2] == 0 || prototypes[i-max2] == 0) ) {
      list1.elem[list1.num_elem] = i;
      list1.num_elem++;
      maps[i] = 0;
    }
    if (list1.num_elem > max_number_in_list) {
      printf("Error list1.num_elem > %d\n",max_number_in_list);
      return 1;
    }
  } 
 
  d = 0;
  while (1) {
  while (list1.num_elem != 0) {
    /* printf("num elem list1: %d\n",list1.num_elem);*/
    for (i=0;i<list1.num_elem;i++) {    
      /* Get element from list1 */
      mapindex = list1.elem[i];
      x = maptox(mapindex,max2);
      y = maptoy(mapindex,max2);      
      new_angle_error = getnewcoordinates(&newx,&newy,x,y,max2,gradientx[y][x],-gradienty[y][x],angle_error);
      newmapindex = newy*max2+newx;
     
      angle_error[newmapindex] = new_angle_error;
      if (newmapindex == mapindex) {
	printf("newmapindex == mapindex!\n");
      }
      if (prototypes[newmapindex] == 2) {
	/* Compute new dstance */
	/* maps[newmapindex] = maps[mapindex] + distance(x,y,newx,newy);*/
	maps[newmapindex] = input_maps[newmapindex];
	if (maps[newmapindex] > 0) {	   	  
	  /* Put new element in list2 */
	  list2.elem[list2.num_elem] = newmapindex;
	  list2.num_elem++;
	  /* Mark the new visited site */
	  prototypes[newmapindex] = 0;
	}	
      }
    }
    /* for (i=0;i<list1.num_elem;i++) {
      mapindex = list1.elem[i];
      for (x=-1;x<2;x++) {  
	for (y=-1;y<2;y++) {
	  if (x==0 && y ==0) continue;  
	  newmapindex = mapindex + x + y*max2;  
	  if (prototypes[newmapindex] == 2) {	      
	    xr = maptox(mapindex,max2);  
	    yr = maptoy(mapindex,max2);    
	    newx = maptox(newmapindex,max2);
	    newy = maptoy(newmapindex,max2);
	    new_angle_error = getnewcoordinates(&provx,&provy,newx,newy,max2,-gradientx[newy][newx],gradienty[newy][newx],angle_error);
            hay que calcular el angulo de error en (newx,newy)
	    provmapindex = provy*max2+provx; 
	    angle_error[newmapindex]=-new_angle_error;
	    if (maps[provmapindex] == -1) {
	      maps[newmapindex] = input_maps[newmapindex];
	    } else {
	      maps[newmapindex] = input_maps[newmapindex];
	    }	  
	   
	    if (maps[newmapindex] != -1) {
	      list3.elem[list3.num_elem] = newmapindex; 
	      list3.num_elem++;	   
	      prototypes[newmapindex] = 0;	  
	      numelemaislados++; 
	    }
	  }
	} 
      } 
    }
    d++;*/
    /* swap(list1.elem,list2.elem); */
    
    index_aux = list1.elem;
    list1.elem = list2.elem;
    list2.elem = index_aux;

    list1.num_elem = list2.num_elem;
    list2.num_elem = 0;
  }
  /* swap(list1.elem,list3.elem); */
  index_aux = list1.elem;
  list1.elem = list3.elem;
  list3.elem = index_aux;
  
  list1.num_elem = list3.num_elem;
  list3.num_elem = 0;
  if (list1.num_elem == 0) break;
  }
  printf("Number of collisions %d\n",counter);
  printf("dmax = %d, numelemaislados %d\n",d,numelemaislados);

  free(angle_error);
  free(list1.elem);
  free(list2.elem);
  free(list3.elem);
  return 0; /* success */
}

/* Mean and standard deviation of the thickness over the band (pixels whose
   input label equals label_cortex). Non-finite map values (e.g. +inf left at a
   boundary pixel) are skipped. Returns the mean; *std gets the standard
   deviation and *npoints the number of band pixels averaged. */
float compute_mean_thickness2D(unsigned char *input, float *maps, int label_cortex, int max1, int max2, int *npoints, float *std) {
  int i, n = 0;
  float mean = 0;

  for (i=0; i<max1*max2; i++) {
    if (input[i] == label_cortex && isfinite(maps[i])) {
      mean += maps[i];
      n++;
    }
  }
  if (n == 0) {
    *npoints = 0;
    *std = 0;
    return 0;
  }
  mean = mean/(float)n;

  *std = 0;
  for (i=0; i<max1*max2; i++) {
    if (input[i] == label_cortex && isfinite(maps[i])) {
      *std += (mean - maps[i])*(mean - maps[i]);
    }
  }
  *std = (n > 1) ? sqrt(*std/(n-1)) : 0;
  *npoints = n;

  return mean;
}
