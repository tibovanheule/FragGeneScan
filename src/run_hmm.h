#ifndef RUN_HMM_H
#define RUN_HMM_H
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "hmm.h"

#include <pthread.h>

typedef struct thread_data {
    FILE *out;
    FILE *aa;
    FILE *dna;
    char *obs_head;
    char *obs_seq;
    int wholegenome;
    int cg;
    int format;
    HMM *hmm;
    TRAIN *train;
} thread_data;

void* thread_func(void *threadarr);

void print_error(const char* error_message, ...);

void print_file_error(const char* error_message, char* file);

void combine(int threadnum,char* out_file,thread_data *threadarr);

#endif