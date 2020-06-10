#ifndef UTIL_LIB_H
#define UTIL_LIB_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

double **dmatrix(int num_row, int num_col);
double *dvector(int nh);
int **imatrix(int num_row, int num_col);
int *ivector(int nh);

void free_dvector(double *v);
void free_dmatrix(double **m,int num_row);
void free_ivector(int *v);
void free_imatrix(int **m,int num_row);

int tr2int (char *nt);
int nt2int (char nt);
int nt2int_rc (char nt);


int trinucleotide (char a, char b, char c);
void get_protein(char *dna, char *protein, int strand, int whole_genome);
void print_usage();
void print_allocation_error(const char *format, ...);

#endif