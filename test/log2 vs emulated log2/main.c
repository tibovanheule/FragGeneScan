#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



double log2_test(double a){
    return log(a)/log(2);
}

double drand(){
    double f=  ((double) rand()*50)/RAND_MAX ;
    return f;
}

//const int test = trinucleotide_pep('G', 'T', 'G');
 #define compile_time_trinucleotide_pep(a,b,c) trinucleotide_pep(a,b,c)
int main(int argc,char *argv[]){
    srand(time(NULL));
    FILE* mathlog, *emulatedlog;
    mathlog = fopen("mathlog.csv","w");
    emulatedlog = fopen("emulatedlog.csv","w");
    for(int i =0; i <= 40000;i++){
        double tester = drand();
        clock_t begin = clock();
        for(int j=0;j<=10000;j++) log2_test(tester);
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        fprintf(emulatedlog,"%f\n",time_spent);

        clock_t begin2 = clock();
        for(int j=0;j<=10000;j++) log2(tester);
        clock_t end2 = clock();
        double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
        fprintf(mathlog,"%f\n",time_spent2);
    }
    fclose(mathlog);
    fclose(emulatedlog);
}

