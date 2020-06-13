#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "util_lib.h"

#ifdef EMULATED_LOG
double log2(double a) {
    return log(a)/log(2);
}
#endif

#define TR_SIZE 14

char* tr_list[TR_SIZE] = { "MM","MI","MD","II","IM","DD","DM","GE","GG","ER","RS","RR","ES","ES1" };

char codon5[5] = { 'A', 'C', 'G', 'T', 'N' };
char codon11[11] = { 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', 'x' };

char codon_code[65] = { 'K','N','K','N',
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
                        'L','F','L','F', 'X'
                      };

char anti_codon_code[65] = { 'F','V','L','I',
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
                             '*','E','Q','K','X'
                           };

/**
* Makes an matrix with datatype double.
* Elements are double pointers en matrix is a double double pointer (**pointer).
* Exits when allocation fails.
*/
double **dmatrix(int num_row, int num_col) {

    double **m = malloc(num_row * sizeof(double*));

    if (!m) print_allocation_error("%s\n", "ERROR: Allocation failure for points to rows in dmatrix()");

    for(int i=0; i<num_row; i++) {
        m[i]=(double *) calloc(num_col, sizeof(double));
        if (!m[i]) print_allocation_error("%s %d %s\n", "ERROR: Allocation failure for the row ", i, " in dmatrix()");
    }
    return m;
}

/**
* Makes an matrix with datatype int.
* Elements are int pointers en matrix is a double int pointer.
* Exits when allocation fails.
*/
int **imatrix(int num_row, int num_col) {

    int i,j;
    int **m = malloc(num_row * sizeof(int*));

    if (!m) print_allocation_error("%s\n", "ERROR: Allocation failure for points to rows in imatrix()");

    for(i=0; i<num_row; i++) {
        m[i]=(int *) calloc(num_col, sizeof(int));
        if (!m[i]) print_allocation_error("%s %d %s\n", "ERROR: Allocation failure for the row ", i," in imatrix()");
    }
    return m;
}

/**
* Makes an vector (array) with datatype double.
* Elements are doubles en vector is a double pointer.
* Exits when allocation fails.
*/
double *dvector(int nh) {

    int j;
    double *v = calloc(nh , sizeof(double));

    if (!v) print_allocation_error("%s\n", "ERROR: Allocation failure in dvector()");
    return v;
}

/**
* Makes an vector array) with datatype int.
* Elements are ints en vector is a int pointer.
* Exits when allocation fails.
*/
int *ivector(int nh) {
    int *v= calloc(nh , sizeof(int));

    if (!v) print_allocation_error("%s\n", "ERROR: Allocation failure in ivector()");

    return v;
}

/**
* Makes an vector array) with datatype int.
* Elements are ints en vector is a int pointer.
* Exits when allocation fails.
*/
int *real_ivector(int* ptr, int nh) {
    int *v= realloc(ptr, nh * sizeof(int));

    if (!v) print_allocation_error("%s\n", "ERROR: Allocation failure in real_ivector()");

    return v;
}

/**
* Frees the memory allocation of an vector with datatype double.
*/
void free_dvector(double *v) {
    free(v);
}

/**
* Frees the memory allocation of an vector with datatype int.
*/
void free_ivector(int *v) {
    free(v);
}

/**
* Frees the memory allocation of an matrix with datatype double.
*/
void free_dmatrix(double **m, int num_row) {
    for(int i=num_row-1; i>=0; i--) free(m[i]);
    free(m);
}

/**
* Frees the memory allocation of an matrix with datatype int.
*/
void free_imatrix(int **m,int num_row) {
    for(int i=num_row-1; i>=0; i--)    free(m[i]);
    free(m);
}


/**
* Converts a given transition to int. Use for example as indexing.
* switch case not possible due the fact that strings are not constonant.
*/
int tr2int(char *tr) {
    size_t index=0;
    while(index < TR_SIZE && (strcmp(tr, tr_list[index])!=0)) ++index;
    return index;
}



int nt2int (char nt) {
    switch(nt) {
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


int nt2int_rc (char nt) {
    switch(nt) {
    case 'A':
    case 'a':
        return 3;
    case 'C':
    case 'c':
        return 2;
    case 'G':
    case 'g':
        return 1;
    case 'T':
    case 't':
        return 0;
    default:
        return 4;
    }
}

int nt2int_rc_indel (char nt) {
    switch(nt) {
    case 'A':
        return 3;
    case 'C':
        return 2;
    case 'G':
        return 1;
    case 'T':
        return 0;
    case 'a':
        return 8;
    case 'c':
        return 7;
    case 'g':
        return 6;
    case 't':
        return 5;
    case 'n':
        return 9;
    case 'x':
        return 10;
    default:
        return 4;
    }
}


int trinucleotide (char a, char b, char c) {

    int freq_id;

    switch(a) {
    case 'A':
    case 'a':
        freq_id = 0;
        break;
    case 'C':
    case 'c':
        freq_id = 16;
        break;
    case 'G':
    case 'g':
        freq_id = 32;
        break;
    case 'T':
    case 't':
        freq_id = 48;
        break;
    default:
        freq_id = 0;
        break;
    }

    switch(b) {
    case 'A':
    case 'a':
        break;
    case 'C':
    case 'c':
        freq_id += 4;
        break;
    case 'G':
    case 'g':
        freq_id += 8;
        break;
    case 'T':
    case 't':
        freq_id += 12;
        break;
    default:
        break;
    }

    switch(c) {
    case 'A':
    case 'a':
        return freq_id;
    case 'C':
    case 'c':
        freq_id += 1;
        return freq_id;
    case 'G':
    case 'g':
        freq_id += 2;
        return freq_id;
    case 'T':
    case 't':
        freq_id += 3;
        return freq_id;
    default:
        return freq_id;
    }
}

int trinucleotide_pep (char a, char b, char c) {

    int freq_id;

    switch(a) {
    case 'A':
    case 'a':
        freq_id = 0;
        break;
    case 'C':
    case 'c':
        freq_id = 16;
        break;
    case 'G':
    case 'g':
        freq_id = 32;
        break;
    case 'T':
    case 't':
        freq_id = 48;
        break;
    default:
        freq_id = 64;
        break;
    }

    if (freq_id < 64) {
        switch(b) {
        case 'A':
        case 'a':
            break;
        case 'C':
        case 'c':
            freq_id += 4;
            break;
        case 'G':
        case 'g':
            freq_id += 8;
            break;
        case 'T':
        case 't':
            freq_id += 12;
            break;
        default:
            freq_id = 64;
            break;
        }
    }

    if (freq_id < 64) {
        switch(c) {
        case 'A':
        case 'a':
            return freq_id;
        case 'C':
        case 'c':
            freq_id += 1;
            return freq_id;
        case 'G':
        case 'g':
            freq_id += 2;
            return freq_id;
        case 'T':
        case 't':
            freq_id += 3;
            return freq_id;
        default:
            return 64;
        }
    }
    return freq_id;
}

/**
* copies dna to dna1 in reverse. and
*/
void get_rc_dna(char *dna, char *dna1) {
    int dna_len = strlen(dna);
    for (int i=0; i<dna_len; i++) dna1[dna_len-i-1] = codon5[nt2int_rc(dna[i])];
}

/**
* copies dna to dna1 in reverse. and
*/
void get_rc_dna_indel(char *dna, char *dna1) {
    int dna_len = strlen(dna);
    for ( int i=0; i<dna_len; i++) dna1[dna_len-i-1] = codon11[nt2int_rc_indel(dna[i])];
}


/**
* Get a protein of dna
* if Whole_genome equals to zero, then we want a short read and stop early.
*/
void get_protein(char *dna, char *protein,  int strand, int whole_genome) {

    int dna_len = strlen(dna);

    if (strand == 1) {
        for (int i=0; i<dna_len; i+=3) protein[i/3] = codon_code[trinucleotide_pep(dna[i], dna[i+1], dna[i+2])];
    } else {
        int protein_len = dna_len/3;

        for (int i=0; i<dna_len; i+=3) {
            protein[(dna_len-i)/3-1] = anti_codon_code[trinucleotide_pep(dna[i], dna[i+1], dna[i+2])];
            protein_len--;
        }
    }

//remove the ending *
    if(protein[strlen(protein) - 1] == '*') protein[strlen(protein) - 1] = 0;

    //alternative start codons still encode for Met
    //E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG and one or two others (e.g., an AUU and possibly a CUG)
    //only consider two major alternative ones, GTG and TTG f

    if(whole_genome == 0) return; //short reads, skip

    if(strand == 1) {
        int s = trinucleotide_pep(dna[0], dna[1], dna[2]);
        if(s == trinucleotide_pep('G', 'T', 'G') || s == trinucleotide_pep('T', 'T', 'G')) protein[0] = 'M';
    } else {
        int s = trinucleotide_pep(dna[dna_len - 3], dna[dna_len - 2], dna[dna_len - 1]);
        if(s == trinucleotide_pep('C', 'A', 'C') || s == trinucleotide_pep('C', 'A', 'A')) protein[0] = 'M';
    }
}

/**
* Print how the program should be used.
* called mainly on help or error.
*/
void print_usage() {

    printf("%s", "USAGE: ./FragGeneScan -s [seq_file_name] -o [output_file_name] -w [1 or 0] -t [train_file_name] (-p [thread_num])\n\n");
    printf("%s", "       Mandatory parameters\n");
    printf("%s", "       [seq_file_name]:    sequence file name including the full path\n");
    printf("%s", "       [output_file_name]: output file name including the full path\n");
    printf("%s", "       [1 or 0]:           1 if the sequence file has complete genomic sequences\n");
    printf("%s", "                           0 if the sequence file has short sequence reads\n");
    printf("%s", "       [train_file_name]:  file name that contains model parameters; this file should be in the \"train\" directory\n");
    printf("%s", "                           Note that four files containing model parameters already exist in the \"train\" directory\n");
    printf("%s", "                           [complete] for complete genomic sequences or short sequence reads without sequencing error\n");
    printf("%s", "                           [sanger_5] for Sanger sequencing reads with about 0.5% error rate\n");
    printf("%s", "                           [sanger_10] for Sanger sequencing reads with about 1% error rate\n");
    printf("%s", "                           [454_5] for 454 pyrosequencing reads with about 0.5% error rate\n");
    printf("%s", "                           [454_10] for 454 pyrosequencing reads with about 1% error rate\n");
    printf("%s", "                           [454_30] for 454 pyrosequencing reads with about 3% error rate\n");
    printf("%s", "                           [illumina_5] for Illumina sequencing reads with about 0.5% error rate\n");
    printf("%s", "                           [illumina_10] for Illumina sequencing reads with about 1% error rate\n\n");
    printf("%s", "       Optional parameter\n");
    printf("%s", "       [thread_num]:       the number of threads used by FragGeneScan; default is 1 thread.\n");
}

/**
* Custom error function to print allocation errors.
* Mostly called from matrix or vector functions.
*/
void print_allocation_error(const char * format,...) {
    va_list list;
    va_start(list,format);
    vfprintf(stderr, format,list);
    va_end(list);
    exit(EXIT_FAILURE);
}