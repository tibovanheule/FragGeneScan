#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

void print_allocation_error(const char * format,...){
    va_list list;
    va_start(list,format);
    vfprintf(stderr, format,list);
    va_end(list);
    exit(EXIT_FAILURE);
}

double *dvector(int nh){
  double *v = (double *) malloc(nh * sizeof(double));

  if (!v) print_allocation_error("%s\n", "ERROR: Allocation failure in dvector()");

  for(int j=0; j<nh; j++) v[j] = 0.0;

  return v;
}


double *dvectorcalloc(int nh){
  double *v = (double *) calloc(nh, sizeof(double));

  if (!v) print_allocation_error("%s\n", "ERROR: Allocation failure in dvector()");
  return v;
}

void free_dvector(double *v){
  free(v);
}

int main(int argc,char *argv[]){

  int num= 35;
  FILE *fp=fopen("malloc.csv", "w");
  for(int k =0;k<20000;k++){
      
      clock_t begin = clock();
      double *m = dvector(num);
      clock_t end = clock();
      if(argc <= 1) for(int i=0; i<num; i++) printf("%s %f", " ", m[i]);
      printf("\n");
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      fprintf(fp,"%f\n",time_spent);
      free(m);
      
  }
  fclose(fp);

  FILE *cp=fopen("calloc.csv", "w");
  for(int k =0;k<20000;k++){

      clock_t begin2 = clock();
      double *c = dvectorcalloc(num);
      clock_t end2 = clock();
  
      if(argc <= 1) for(int j=0; j<num; j++) printf("%s %f", " ", c[j]);
      printf("\n");
      
      double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
      fprintf(cp,"%f\n",time_spent2);

      free(c);
  }
  fclose(cp);
}

