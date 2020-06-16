#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



int tr2int (char nt){

  int result;

  if      (nt == 'A' || nt == 'a'){  result = 0; }
  else if (nt == 'C' || nt == 'c'){  result = 1; }
  else if (nt == 'G' || nt == 'g'){  result = 2; }
  else if (nt == 'T' || nt == 't'){  result = 3; }
  else                            {  result = 4; }

  return result;
}

int tr2ints (char nt){

  switch(nt){
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
        return 3;
    default:
        return 4;
  }
}

int main(int argc,char *argv[]){

  FILE *fp = fopen("iffen.csv", "w");
  for(int k =0;k<20000;k++){
      
      clock_t begin = clock();
        tr2int('A');
        tr2int('a');
        tr2int('c');
        tr2int('g');
        tr2int('t');
        tr2int('T');
        tr2int('G');

        tr2int('C');
tr2int('h');
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      fprintf(fp,"%f\n",time_spent);
      
  }
  fclose(fp);

  FILE *cp=fopen("switch.csv", "w");
  for(int k =0;k<20000;k++){

      clock_t begin2 = clock();
        tr2ints('A');
        tr2ints('a');
        tr2ints('c');
        tr2ints('g');
        tr2ints('t');
        tr2ints('T');
        tr2ints('G');

        tr2ints('C');
tr2ints('h');
      clock_t end2 = clock();
      double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
      fprintf(cp,"%f\n",time_spent2);
        printf("%d\n",tr2int('A')==tr2ints('a'));
        printf("%d\n",tr2int('a')==tr2ints('A'));
        printf("%d\n",tr2int('T')==tr2ints('t'));
        printf("%d\n",tr2int('t')==tr2ints('T'));
        printf("%d\n",tr2int('G')==tr2ints('G'));
        printf("%d\n",tr2int('g')==tr2ints('g'));
        printf("%d\n",tr2int('C')==tr2ints('C'));

        printf("%d\n",tr2int('c')==tr2ints('c'));
printf("%d\n",tr2int('h')==tr2ints('h'));

  }
  fclose(cp);
}

