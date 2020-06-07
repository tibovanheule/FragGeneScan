#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



/**
* Makes an matrix with datatype double.
* Elements are double pointers en matrix is a double double pointer (**pointer).
* Exits when allocation fails.
*/
double **dmatrix(int num_row, int num_col){

  int i, j;
  double **m = (double **) malloc(num_row * sizeof(double*));

  if (!m) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in dmatrix()");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<num_row; i++) {
    m[i]=(double *) malloc(num_col * sizeof(double));
    if (!m[i]) {
      fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i, " in dmatrix()");
      exit(EXIT_FAILURE);
    }

    for(j=0; j<num_col; j++) m[i][j] = 0.0;
    
  }
  return m;
}


double **dmatrixcalloc(int num_row, int num_col){

  int i, j;
  double **m = (double **) malloc(num_row * sizeof(double*));

  if (!m) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in dmatrix()");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<num_row; i++) {
    m[i]=(double *) calloc(num_col,sizeof(double));
    if (!m[i]) {
      fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i, " in dmatrix()");
      exit(EXIT_FAILURE);
    }    
  }
  return m;
}

void free_dmatrix(double **m, int num_row){
  for(int i=num_row-1; i>=0; i--) free(m[i]);
  free(m);
}

int main(int argc,char *argv[]){

  int num_row= 4;
  int num_col = 4;
  FILE *fp=fopen("malloc.csv", "w");
  for(int k =0;k<20000;k++){
      
      clock_t begin = clock();
      double **m = dmatrix(num_row, num_col);
      clock_t end = clock();
      if(argc <= 1) for(int i=0; i<num_row; i++) {
        for(int j=0; j<num_col; j++) printf("%s %f", " ", m[i][j]);
        printf("\n");
      }
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      fprintf(fp,"%f\n",time_spent);
      free_dmatrix(m,num_row);
      
  }
  fclose(fp);

  FILE *cp=fopen("calloc.csv", "w");
  for(int k =0;k<20000;k++){

      clock_t begin2 = clock();
      double **c = dmatrixcalloc(num_row, num_col);
      clock_t end2 = clock();
  
      if(argc <= 1) for(int i=0; i<num_row; i++) {
        for(int j=0; j<num_col; j++) printf("%s %f", " ", c[i][j]);
        printf("\n");
      }

      double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
      fprintf(cp,"%f\n",time_spent2);

      free_dmatrix(c,num_row);
  }
  fclose(cp);
}

