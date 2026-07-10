/* Copyright (c) Ruben Cardenes Almeida 22/03/2002 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "thickness3D.h"
#include "laplace3D.h"

#define PI 3.1415927
#define INF 9999999
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

float distance(int x1,int y1,int x2,int y2) {
  return sqrt ((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
}

float distanceYezzi_reverse3D_relax(float ***gradientx,float ***gradienty,float ***gradientz,int newmapindex, int x, int y,int z,float *maps,int max1,int max2, float r,float hx,float hy,float hz) {
  float distf;
  int flag = 0;
  
  distf = 1.0;
  if (gradientx[z][y][x] > 0) { 
    if (maps[newmapindex+1] > -1) {
      distf += fabs(gradientx[z][y][x])*maps[newmapindex+1]/hx;
      flag = 1;
    }
  } else if(gradientx[z][y][x] < 0) {
    if (maps[newmapindex-1] > -1) {
      distf += fabs(gradientx[z][y][x])*maps[newmapindex-1]/hx;
      flag = 1;
    }
  }

  if (gradienty[z][y][x] > 0) {
    if (maps[newmapindex+max2] > -1) { 
      distf += fabs(gradienty[z][y][x])*maps[newmapindex+max2]/hy;
      flag = 1;
    }
  } else if(gradienty[z][y][x] < 0)  {
    if (maps[newmapindex-max2] > -1) { 
      distf += fabs(gradienty[z][y][x])*maps[newmapindex-max2]/hy;
      flag = 1;
    }    
  }

  if (gradientz[z][y][x] > 0) {
    if (maps[newmapindex+max1*max2] > -1) {
      distf += fabs(gradientz[z][y][x])*maps[newmapindex+max1*max2]/hz;
      flag = 1;
    }
  } else if (gradientz[z][y][x] < 0) {
    if (maps[newmapindex-max1*max2] > -1) {
      distf += fabs(gradientz[z][y][x])*maps[newmapindex-max1*max2]/hz;
      flag = 1;
    }
  }

  if (flag == 0) {
    distf = -1;
  } else {
    distf = distf/(fabs(gradientx[z][y][x])/hx + fabs(gradienty[z][y][x])/hy + fabs(gradientz[z][y][x])/hz);   
  }

  return distf;
}

float distanceYezzi3D_relax(float ***gradientx,float ***gradienty,float ***gradientz,int newmapindex,int x, int y,int z, float *maps,int max1,int max2,float r,float hx, float hy, float hz) {
  float distf;
  int flag = 0;

  distf = 1.0;
  if (gradientx[z][y][x] < 0) { 	      
    if (maps[newmapindex+1] > -1) {
      flag = 1;
      distf += fabs(gradientx[z][y][x])*maps[newmapindex+1]/hx;
    }    	     
  } else {
    if (maps[newmapindex-1] > -1) {     
      flag = 1;
      distf += fabs(gradientx[z][y][x])*maps[newmapindex-1]/hx;
    }
  }
  if (gradienty[z][y][x] < 0) {    
    if (maps[newmapindex+max2] > -1) {
      distf += fabs(gradienty[z][y][x])*maps[newmapindex+max2]/hy;
      flag = 1;
    }
  } else {    
    if (maps[newmapindex-max2] > -1) {
      distf += fabs(gradienty[z][y][x])*maps[newmapindex-max2]/hy;
      flag = 1;
    }
  }

  if (gradientz[z][y][x] < 0) {
    if (maps[newmapindex+max1*max2] > -1) {
      distf += fabs(gradientz[z][y][x])*maps[newmapindex+max1*max2]/hz;
      flag = 1;
    }
  } else {
    if (maps[newmapindex-max1*max2] > -1) {
      distf += fabs(gradientz[z][y][x])*maps[newmapindex-max1*max2]/hz;
      flag = 1;
    }
  }

  if (flag == 0) {
    distf = -1;
  } else {
    distf = distf/(fabs(gradientx[z][y][x])/hx + fabs(gradienty[z][y][x])/hy + fabs(gradientz[z][y][x])/hz);
  }
  return distf;
}

float distanceYezzi_reverse3D(float ***gradientx,float ***gradienty,float ***gradientz,int newmapindex, int x, int y,int z,float *maps,int max1,int max2, float r,float hx,float hy,float hz) {
  float distf;

  if (gradientx[z][y][x] > 0) { 
    if (maps[newmapindex+1] == -1 && fabs(gradientx[z][y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex+1] == -1) { 
	distf = 1.0;
      } else {
	distf = 1.0 + fabs(gradientx[z][y][x])*maps[newmapindex+1]/hx;	      
      }
    }
  } else {
    if (maps[newmapindex-1] == -1 && fabs(gradientx[z][y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex-1] == -1) { 
	distf = 1.0;
      } else {
	distf = 1.0 + fabs(gradientx[z][y][x])*maps[newmapindex-1]/hx;
      }
    }
  }
  if (gradienty[z][y][x] > 0) {
    if (maps[newmapindex+max2] == -1 && fabs(gradienty[z][y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex+max2] > -1) { 
	distf += fabs(gradienty[z][y][x])*maps[newmapindex+max2]/hy;
      }
    }
  } else {
    if (maps[newmapindex-max2] == -1 && fabs(gradienty[z][y][x]) > r ) {
      return -1;
    } else {
      if (maps[newmapindex-max2] > -1) { 
	distf += fabs(gradienty[z][y][x])*maps[newmapindex-max2]/hy;
      }
    }
  }

  if (gradientz[z][y][x] > 0) {
    if (maps[newmapindex+max1*max2] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex+max1*max2] > -1) {
	distf += fabs(gradientz[z][y][x])*maps[newmapindex+max1*max2]/hz;
      }
    }
  } else {
    if (maps[newmapindex-max1*max2] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    }  else {
      if (maps[newmapindex-max1*max2] > -1 ) {
	distf += fabs(gradientz[z][y][x])*maps[newmapindex-max1*max2]/hz;
      }
    }
  }
 
  distf = distf/(fabs(gradientx[z][y][x])/hx + fabs(gradienty[z][y][x])/hy + fabs(gradientz[z][y][x])/hz);
 
  return distf;
}

float distanceYezzi3D(float ***gradientx,float ***gradienty,float ***gradientz,int newmapindex,int x, int y,int z, float *maps,int max1,int max2,float r,float hx, float hy, float hz) {
  float distf;
  if (gradientx[z][y][x] < 0) { 	      
    if (maps[newmapindex+1] == -1 && fabs(gradientx[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex+1] == -1) {
	distf = 1.0;
      } else {
	distf = 1.0 + fabs(gradientx[z][y][x])*maps[newmapindex+1]/hx;
      }
    }	      
  } else {
    if (maps[newmapindex-1] == -1 && fabs(gradientx[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex-1] == -1) {
	distf = 1.0;
      } else { 
	distf = 1.0 + fabs(gradientx[z][y][x])*maps[newmapindex-1]/hx;
      }
    }
  }
  if (gradienty[z][y][x] < 0) {
    if (maps[newmapindex+max2] == -1 && fabs(gradienty[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex+max2] > -1) {
	distf += fabs(gradienty[z][y][x])*maps[newmapindex+max2]/hy;
      }
    }
  } else {
    if (maps[newmapindex-max2] == -1 && fabs(gradienty[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex-max2] > -1) {
	distf += fabs(gradienty[z][y][x])*maps[newmapindex-max2]/hy;
      }
    }
  }

  if (gradientz[z][y][x] < 0) {
    if (maps[newmapindex+max1*max2] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex+max1*max2] > -1) {
	distf += fabs(gradientz[z][y][x])*maps[newmapindex+max1*max2]/hz;
      }
    }
  } else {
    if (maps[newmapindex-max1*max2] == -1 && fabs(gradientz[z][y][x]) > r) {
      return -1;
    } else {
      if (maps[newmapindex-max1*max2] > -1) {
	distf += fabs(gradientz[z][y][x])*maps[newmapindex-max1*max2]/hz;
      }
    }
  }

  distf = distf/(fabs(gradientx[z][y][x])/hx + fabs(gradienty[z][y][x])/hy + fabs(gradientz[z][y][x])/hz);
  return distf;
}

int thickness3DYezzi(unsigned char* prototypes,int max1, int max2, int max3, float *maps, float*** laplacefield,float*** gradientx, float*** gradienty, float*** gradientz, int num_it, float hx, float hy, float hz) {
  int i,j,x,y,z,xr,yr,zr,mapindex,newmapindex,aux,d,l,flag;
  int *index_aux;
  float distf,r;
  struct list {
    int num_elem;
    int *elem;
  } list1; 
  struct list list2,list_aux;
  int max_number_in_list = 300000;
  unsigned char* prot_copia;
  
  prot_copia = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2*max3);

  for (j=0;j<max1*max2*max3;j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list1.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list1.num_elem = 0;

  list2.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list2.num_elem = 0;
  
  list_aux.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list_aux.num_elem = 0;
  printf("iteration ");
  for (l=0;l<num_it;l++) {
    d = 0;
    if (l == 0) {
      r=0.3;
    } else {
      r=0.0;
    }
    printf("%02d\b\b",l);
    fflush(0);
    for (i=0;i<max1*max2*max3;i++) {
      prot_copia[i] =prototypes[i]; /* to recover the domain after the first iteration */
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
	for (z=-1;z<2;z++) {
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      if ((abs(x) + abs(y) + abs(z))!=1) continue; /* 3D neighborhood of size 6*/
	      newmapindex = mapindex + max1*max2*z + y*max2 + x;
	      xr = maptox3d(newmapindex,max1,max2);
	      yr = maptoy3d(newmapindex,max1,max2);
	      zr = maptoz3d(newmapindex,max1,max2);
	      if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1 || zr < 0 || zr >= max3) continue;
	      if (prot_copia[newmapindex] == 2) {
		/* Compute new distance */
		distf = distanceYezzi3D(gradientx,gradienty,gradientz,newmapindex,
					xr,yr,zr,maps,max1,max2,r,hx,hy,hz);

		if (distf > 0 && distf < INF) {
		  if (fabs(maps[mapindex] - distf) < 1) {
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
	for (z=-1;z<2;z++) {
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      if ((abs(x) + abs(y) + abs(z))!=1) continue; /* 3D neighborhood of size 6*/
	      newmapindex = mapindex + max1*max2*z + y*max2 + x;
	      xr = maptox3d(newmapindex,max1,max2);
	      yr = maptoy3d(newmapindex,max1,max2);
	      zr = maptoz3d(newmapindex,max1,max2);
	      if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1 || zr < 0 || zr >= max3) continue;
	      if (prot_copia[newmapindex] == 2) {
		/* Compute new distance */
		if (l==0) {
		  distf = distanceYezzi3D_relax(gradientx,gradienty,gradientz,newmapindex,
						xr,yr,zr,maps,max1,max2,r,hx,hy,hz);
		} else {
		  distf = distanceYezzi3D(gradientx,gradienty,gradientz,newmapindex,
					  xr,yr,zr,maps,max1,max2,r,hx,hy,hz);
		}
		if (distf > 0 && distf < INF) {
		  if (l==0 || fabs(maps[mapindex] - distf) < 1) {
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
  printf("\nmaximum bucket = %d\n",d);

  free(list1.elem);
  free(list2.elem);
  free(list_aux.elem);
  free(prot_copia);
  return 0; /* success */
}

int thickness3DYezzi_reverse(unsigned char* prototypes,int max1, int max2, int max3, float *maps, float*** laplacefield,float*** gradientx, float*** gradienty, float*** gradientz, int num_it, float hx, float hy, float hz) {
  int i,j,x,y,z,xr,yr,zr,mapindex,newmapindex,chekmapindex,reject,aux,d,flag,l,ag;
  int *index_aux;
  float distf,r;
  struct list {
    int num_elem;
    int *elem;
  } list1; 
  struct list list2,list_aux;
  int max_number_in_list = 300000;
  unsigned char* prot_copia;
  FILE *fp;

  prot_copia = (unsigned char*)malloc(sizeof(unsigned char)*max1*max2*max3);
  
  for (j=0;j<max1*max2*max3;j++) {
    maps[j] = -1;
    prot_copia[j] = prototypes[j];
  }

  list1.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list1.num_elem = 0;
  
  list2.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list2.num_elem = 0;
  
  list_aux.elem = (int*)malloc(sizeof(int)*max_number_in_list);
  list_aux.num_elem = 0;
  printf("iteration ");
  for (l=0;l<num_it;l++) {
    d = 0;
    if (l == 0) {
      r=0.3;
    } else {
      r=0.0;
    }
    printf("%02d\b\b",l);
    fflush(0);
    for (i=0;i<max1*max2*max3;i++) {
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
	for (z=-1;z<2;z++) {
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      if ((abs(x) + abs(y) + abs(z))!=1) continue; /* 3D neighborhood of size 6*/
	      newmapindex = mapindex + max1*max2*z + y*max2 +x;
	      xr = maptox3d(newmapindex,max1,max2);
	      yr = maptoy3d(newmapindex,max1,max2);
	      zr = maptoz3d(newmapindex,max1,max2);
	      if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1 || zr < 0 || zr >= max3) continue;
	      if (prot_copia[newmapindex] == 2) {
		/* Compute new distance */
		distf = distanceYezzi_reverse3D(gradientx,gradienty,gradientz,newmapindex,
						xr,yr,zr,maps,max1,max2,r,hx,hy,hz);

		if (distf > 0 && distf < INF) {
		  if (fabs(maps[mapindex] - distf) < 1.5) {
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
	for (z=-1;z<2;z++) {
	  for (x=-1;x<2;x++) {
	    for (y=-1;y<2;y++) {
	      if ((abs(x) + abs(y) + abs(z))!=1) continue; /* 3D neighborhood of size 6*/
	      newmapindex = mapindex + max1*max2*z + y*max2 + x;
	      xr = maptox3d(newmapindex,max1,max2);
	      yr = maptoy3d(newmapindex,max1,max2);
	      zr = maptoz3d(newmapindex,max1,max2);
	      if (xr < 0 || xr >= max2 || yr < 0 || yr >= max1 || zr < 0 || zr >= max3) continue;
	      /*if (prototypes[newmapindex] == 2 || prototypes[newmapindex] == 0) {*/
	      if (prot_copia[newmapindex] == 2) {
		/* Compute new distance */
		if (l==0) {
		  distf = distanceYezzi_reverse3D_relax(gradientx,gradienty,gradientz,newmapindex,
							xr,yr,zr,maps,max1,max2,r,hx,hy,hz);
		} else {
		  distf = distanceYezzi_reverse3D(gradientx,gradienty,gradientz,newmapindex,
						  xr,yr,zr,maps,max1,max2,r,hx,hy,hz);
		}
		if (distf > 0 && distf < INF) {
		  if (fabs(maps[mapindex] - distf) < 1.5 || l ==0) {
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
  printf("\nmaximum bucket = %d\n",d);
 
  /* codigo control */
  d = 0;ag=0;
  for (i=0;i<max1*max2*max3;i++) {
    if (maps[i] > num_it ) {
      d++;
    }
    if (prototypes[i] == 2 && maps[i] == -1) {
      ag++;
    }
  }
  /* printf("Num potenciales errores %d, agujeros %d\n",d,ag);*/
  /* */

  free(list1.elem);
  free(list2.elem);
  free(list_aux.elem);
  free(prot_copia);
  return 0; /* success */
}

float compute_mean_thickness(unsigned char *input,float *maps, int label_cortex,int boundary_l, int max1, int max2, int max3, float *std ) {
  int i,j,k,sum = 0, npoints = 0;
  float mean = 0;

  for(k=0; k<max3; k++) {
    for(j=0; j<max2; j++) {
      for(i=0; i<max1; i++) {     
	if ((i==0)||(j==0)||(k==0)||(i==max1-1)||(j==max2-1)||(k==max3-1) ) {
	  /* nothing to do */	 
	} else if ( (input[sum]==label_cortex) &&
		    ((input[sum+1]==boundary_l)||
		     (input[sum-1]==boundary_l)||
		     
		     (input[sum+max1]==boundary_l)||
		     (input[sum-max1]==boundary_l)||
		     
		     (input[sum+max1*max2]==boundary_l)||
		     (input[sum-max1*max2]==boundary_l))) {	  
	  mean += maps[sum];
	  npoints++;	  
	}
	sum++;
      }
    }
  }
  
  mean = mean/(float)npoints;
  sum=0;
  (*std) = 0;
  for(k=0; k<max3; k++) {
    for(j=0; j<max2; j++) {
      for(i=0; i<max1; i++) {     
	if ((i==0)||(j==0)||(k==0)||(i==max1-1)||(j==max2-1)||(k==max3-1) ) {
	  /* nothing to do */	 
	} else if ( (input[sum]==label_cortex) &&
		    ((input[sum+1]==boundary_l)||
		     (input[sum-1]==boundary_l)||
		     
		     (input[sum+max1]==boundary_l)||
		     (input[sum-max1]==boundary_l)||
		     
		     (input[sum+max1*max2]==boundary_l)||
		     (input[sum-max1*max2]==boundary_l))) {	  
	  (*std) += (mean - maps[sum])*(mean - maps[sum]);
	}
	sum++;
      }
    }
  }

  printf("npoints %d\n",npoints);
  (*std) = sqrt((*std)/(npoints -1));

  return mean;
}

float compute_mean_thickness_volume(unsigned char *input,float *maps, int label_cortex, int max1, int max2, int max3,float *std) {
  int i, npoints = 0;
  float mean = 0;

  for(i=0; i<max1*max2*max3; i++) {     
    if (input[i] == label_cortex) {	  
      mean += maps[i];
      npoints++;
    }	
  }
  mean = mean/(float)npoints;

  (*std) = 0;
  for(i=0; i<max1*max2*max3; i++) {     
    if (input[i] == label_cortex) {	  
      (*std) += (mean - maps[i])*(mean - maps[i]);
    }
  }

  printf("npoints %d\n",npoints);
  (*std) = sqrt((*std)/(npoints -1));

  return mean;
}
