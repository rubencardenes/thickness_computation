#include <stdio.h>

int im2vtk(unsigned char *data,int height,int width, int depth, char* fichero) {
  int z,x,y,i;
  FILE * fd;
  
  /* Prepare header for vtk file */
   fd = fopen(fichero,"w");
   fprintf(fd,"# vtk DataFile Version 2.0\n");
   fprintf(fd,"Skeleton\n");
   fprintf(fd,"ASCII\n");
   fprintf(fd,"DATASET STRUCTURED_POINTS\n");
   fprintf(fd,"DIMENSIONS %d %d %d\n", height, width, depth);
   fprintf(fd,"ASPECT_RATIO 1 1 1\n");
   fprintf(fd,"ORIGIN 0 0 0\n");   
   fprintf(fd,"POINT_DATA %d\n", height*width*depth);
   fprintf(fd,"SCALARS volume_scalars float\n");
   fprintf(fd,"LOOKUP_TABLE default\n");

   /* write the image data in ascii vtk*/
   i = 0;
   for (z=0; z<depth; z++) {
     for(x=0; x<height; x++) {
       for(y=0; y<width; y++) {
	  fprintf(fd,"%f\n",(float)data[i]);
	  i++;
        } 
     }
   }
   fclose(fd);
   printf("converting .im to .vtk finished\n");
   return 0;
}

int main(int argc, char* argv[]) {
  unsigned short* data;
  FILE* fp;
  char fichero[200];
  int i,j,k,nrows,ncols,nslices;
  float radius = 10;
  float rx = 15, ry = 35, rz = 20;
  float xc = 40,yc = 40,zc = 40; 
  float el,c;

  nrows = 80;
  ncols = 80;
  nslices = 80;
  sprintf(fichero,"%s","phantom_sphere.vols");
  
  data=(unsigned short*)malloc(sizeof(unsigned short)*nrows*ncols*nslices);
  for (k=0;k<nslices;k++) {
    for (j=0;j<nrows;j++) {
      for (i=0;i<ncols;i++) {
	el = (float)((i-xc)*(i-xc))/(rx*rx) + ((j-yc)*(j-yc))/(ry*ry) + ((k-zc)*(k-zc))/(rz*rz);
	if (el > 1) {
	  data[nrows*ncols*k+ncols*j+i] = 255;
	} else {
	  c =  (float)(i-xc)*(i-xc) + (j-yc)*(j-yc) + (k-zc)*(k-zc);
	  if ( c > radius*radius ) {
	    data[nrows*ncols*k+ncols*j+i] = 2;
	  } else {	 	  
	    data[nrows*ncols*k+ncols*j+i] = 3;
	  }
	}
      }
    }
  }
  fp=fopen(fichero,"w");
  fwrite(data,sizeof(unsigned short),ncols*nrows*nslices,fp);
  fclose(fp);

  free(data);
  return 0;
}
