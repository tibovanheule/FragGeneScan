#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



int tr2int (char *tr){

  int result;

  if      (strcmp(tr, "MM")==0){   result = 0; }
  else if (strcmp(tr, "MI")==0){   result = 1; }
  else if (strcmp(tr, "MD")==0){   result = 2; }
  else if (strcmp(tr, "II")==0){   result = 3; }
  else if (strcmp(tr, "IM")==0){   result = 4; }
  else if (strcmp(tr, "DD")==0){   result = 5; }
  else if (strcmp(tr, "DM")==0){   result = 6; }
  else if (strcmp(tr, "GE")==0){   result = 7; }
  else if (strcmp(tr, "GG")==0){   result = 8; }
  else if (strcmp(tr, "ER")==0){   result = 9; }
  else if (strcmp(tr, "RS")==0){   result = 10;}
  else if (strcmp(tr, "RR")==0){   result = 11;}
  else if (strcmp(tr, "ES")==0){   result = 12;}    /* ES: E+ -> S+, E- -> S- */
  else if (strcmp(tr, "ES1")==0){   result = 13;}   /* ES1: E+ -> S-, E- -> S+ */

  return result;
}
#define TR_SIZE 14
int tr2ints(char *tr){

  size_t index=0;
  char* list[TR_SIZE] = {"MM","MI","MD","II","IM","DD","DM","GE","GG","ER","RS","RR","ES","ES1"};
  while(index < TR_SIZE && list[index] != tr) ++index;
  return index;
}

int main(int argc,char *argv[]){

  FILE *fp = fopen("iffen.csv", "w");
  for(int k =0;k<20000;k++){
      
      clock_t begin = clock();
        tr2int("MM");
        tr2int("MI");
        tr2int("MD");
        tr2int("II");
        tr2int("IM");
        tr2int("DD");
        tr2int("DM");

        tr2int("GE");
        tr2int("GG");
        tr2int("ER");
        tr2int("RS");
        tr2int("RR");
        tr2int("ES");
        tr2int("ES1"); 
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      fprintf(fp,"%f\n",time_spent);
      
  }
  fclose(fp);

  FILE *cp=fopen("array.csv", "w");
  for(int k =0;k<20000;k++){

      clock_t begin2 = clock();
        tr2ints("MM");
        tr2ints("MI");
        tr2ints("MD");
        tr2ints("II");
        tr2ints("IM");
        tr2ints("DD");
        tr2ints("DM");

        tr2ints("GE");
        tr2ints("GG");
        tr2ints("ER");
        tr2ints("RS");
        tr2ints("RR");
        tr2ints("ES");
        tr2ints("ES1"); 
      clock_t end2 = clock();
int test = tr2ints("ES1");
      double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
      fprintf(cp,"%f\n",time_spent2);
        printf("%d\n",tr2int("MM")==tr2ints("MM"));
        printf("%d\n",tr2int("MI")==tr2ints("MI"));
        printf("%d\n",tr2int("MD")==tr2ints("MD"));
        printf("%d\n",tr2int("II")==tr2ints("II"));
        printf("%d\n",tr2int("IM")==tr2ints("IM"));
        printf("%d\n",tr2int("DD")==tr2ints("DD"));
        printf("%d\n",tr2int("DM")==tr2ints("DM"));

        printf("%d\n",tr2int("GE")==tr2ints("GE"));
        printf("%d\n",tr2int("GG")==tr2ints("GG"));
        printf("%d\n",tr2int("ER")==tr2ints("ER"));
        printf("%d\n",tr2int("RS")==tr2ints("RS"));
        printf("%d\n",tr2int("RR")==tr2ints("RR"));
        printf("%d\n",tr2int("ES")==tr2ints("ES"));
        printf("%d\n",tr2int("ES1")==tr2ints("ES1")); 

  }
  fclose(cp);
}

