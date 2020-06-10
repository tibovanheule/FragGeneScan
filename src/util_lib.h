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

#define TR_SIZE 14
const char* tr_list[TR_SIZE] = {"MM","MI","MD","II","IM","DD","DM","GE","GG","ER","RS","RR","ES","ES1"};

const char codon5[5] = {'A', 'C', 'G', 'T', 'N'};
const char codon11[11] = {'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', 'x'};

const char codon_code[65] = {'K','N','K','N',
			 'T','T','T','T',
			 'R','S','R','S',
			 'I','I','M','I',
			 'Q','H','Q','H',
			 'P','P','P','P',
			 'R','R','R','R',
			 'L','L','L','L',
			 'E','D','E','D',
			 'A','A','A','A',
			 'G','G','G','G',
			 'V','V','V','V',
			 '*','Y','*','Y',
			 'S','S','S','S',
			 '*','C','W','C',
			 'L','F','L','F', 'X'};

const char anti_codon_code[65] = {'F','V','L','I',
			 'C','G','R','S',
			 'S','A','P','T',
			 'Y','D','H','N',
			 'L','V','L','M',
			 'W','G','R','R',
			 'S','A','P','T',
			 '*','E','Q','K',
			 'F','V','L','I',
			 'C','G','R','S',
			 'S','A','P','T',
			 'Y','D','H','N',
			 'L','V','L','I',
			 '*','G','R','R',
			 'S','A','P','T',
			 '*','E','Q','K','X'};

#endif