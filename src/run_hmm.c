#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "run_hmm.h"
#include "util_lib.h"

#include "thpool.h"

#define ADD_LEN 1024
#define STRINGLEN 4096+1
#define DELIMI " "
/**
* Entry point of program
* 1. Initialization of variables and datatypes
* 2. Check File acessiblity
*/
int main (int argc, char **argv) {
    /*Initialization datastructures*/
    TRAIN *train = malloc(sizeof(TRAIN));
    HMM *hmm = malloc(sizeof(HMM));

    clock_t start = clock();

    /* number of threads */
    int threadnum = 1;

    char hmm_file[STRINGLEN] = "";
    char out_header[STRINGLEN] = "";

    char seq_file[STRINGLEN] = "";

    char train_file[STRINGLEN] = "";
    char mstate_file[STRINGLEN] = "";
    char rstate_file[STRINGLEN] = "";
    char nstate_file[STRINGLEN] = "";
    char sstate_file[STRINGLEN] = "";
    char pstate_file[STRINGLEN] = "";
    char s1state_file[STRINGLEN] = "";     /* stop codon of gene in - stand */
    char p1state_file[STRINGLEN] = "";
    char dstate_file[STRINGLEN] = "";
    char train_dir[STRINGLEN] = "";
    char mystring[STRINGLEN] = "";




    getcwd(train_dir,STRINGLEN);
    strcat(train_dir, "/train/");
    strcpy(mstate_file, train_dir);
    strcat(mstate_file, "gene");
    strcpy(rstate_file, train_dir);
    strcat(rstate_file, "rgene");
    strcpy(nstate_file, train_dir);
    strcat(nstate_file, "noncoding");
    strcpy(sstate_file, train_dir);
    strcat(sstate_file, "start");
    strcpy(pstate_file, train_dir);
    strcat(pstate_file, "stop");
    strcpy(s1state_file, train_dir);
    strcat(s1state_file, "stop1");
    strcpy(p1state_file, train_dir);
    strcat(p1state_file, "start1");
    strcpy(dstate_file, train_dir);
    strcat(dstate_file, "pwm");


    /* If there is less then 9 arguments, then we halt because not everything is given. */
    if (argc <= 8) {
        print_error("ERROR: You missed some parameters for input\n");
    }

    // arguments variables
    int wholegenome, format=0;

    // variables
    int c;
    while ((c=getopt(argc, argv, "fs:o:w:t:p:")) != -1) {
        switch (c) {
        case 's':
            strcpy(seq_file, optarg);
            if (access(seq_file, F_OK)==-1) print_error("ERROR: Sequence file [%s] does not exist\n",seq_file);
            break;
        case 'w':
            wholegenome = atoi(optarg);
            if (wholegenome != 0 && wholegenome != 1)	print_error("ERROR: An incorrect value for the option -w was entered\n");
            break;
        case 'p':
            threadnum = atoi(optarg);
            if (threadnum < 1) print_error("ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
            printf("Using %d threads.\n", threadnum);
            break;
        case 'o':
            strcpy(out_header, optarg);
            break;
        case 't':
            strcpy(train_file, optarg);
            strcpy(hmm_file, train_dir);
            strcat(hmm_file, train_file);
            break;
        case 'f':
            format = 1;
            break;
        }
    }


    /* check whether the specified files exist */
    if (access(mstate_file, F_OK)==-1) print_file_error("Forward prob. file [%s] does not exist\n", mstate_file);
    if (access(rstate_file, F_OK)==-1) print_file_error("Backward prob. file [%s] does not exist\n", rstate_file);
    if (access(nstate_file, F_OK)==-1) print_file_error("noncoding prob. file [%s] does not exist\n", nstate_file);
    if (access(sstate_file, F_OK)==-1) print_file_error("start prob. file [%s] does not exist\n", sstate_file);
    if (access(pstate_file, F_OK)==-1) print_file_error("stop prob. file [%s] does not exist\n", pstate_file);
    if (access(s1state_file, F_OK)==-1) print_file_error("start1 prob. file [%s] does not exist\n", s1state_file);
    if (access(p1state_file, F_OK)==-1) print_file_error("stop1 prob. file [%s] does not exist\n", p1state_file);
    if (access(dstate_file, F_OK)==-1) print_file_error("pwm dist. file [%s] does not exist\n", dstate_file);
    if (access(hmm_file, F_OK)==-1) print_file_error("hmm file [%s] does not exist\n", hmm_file);

    /* read all initial model */
    /* Read the given files and give HMM and TRAIN */
    get_train_from_file(hmm_file, hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, train);

    // Initialize thread data structure
    thread_data *threadarr = malloc(sizeof(thread_data)* threadnum);

    for (int i = 0; i < threadnum; i++)   {
        if(threadnum > 1) {
            sprintf(mystring, "%s.out.tmp.%d", out_header, i);
            threadarr[i].out = fopen(mystring, "w");
            sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
            threadarr[i].aa = fopen(mystring, "w");
            sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
            threadarr[i].dna = fopen(mystring, "w");
        } else {
            sprintf(mystring, "%s.out", out_header);
            threadarr[i].out = fopen(mystring, "w");
            sprintf(mystring, "%s.faa", out_header);
            threadarr[i].aa = fopen(mystring, "w");
            sprintf(mystring, "%s.ffn", out_header);
            threadarr[i].dna = fopen(mystring, "w");
        }

        threadarr[i].hmm = hmm;
        threadarr[i].train = train;

        threadarr[i].wholegenome = wholegenome;
        threadarr[i].format = format;
    }


    // used for th sum of bpcounts
    int sequence_offset = 0;
    // totaal
    int  total=0;
    // array of threads (variables needed for pthreads)
    pthread_t thread[threadnum];
    //open file handler
    FILE *fp = fopen (seq_file, "r");
    // initial size of sequentie string, this will increase if there are strings in the
    size_t size = 150;
    char* sequentie = malloc(size* sizeof(char));

    threadpool thpool = thpool_init(threadnum);
    while ( fgets (mystring, sizeof mystring, fp)  ) {
        if (mystring[0] == '>' || feof(fp)) {

            if (total > 0) {
                threadarr[(total-1) % threadnum].obs_seq = strdup(sequentie);
                sequence_offset = 0;
            }
            if (total > 0 && (total % threadnum == 0 || feof(fp))) {
                // Deal with the thread
                for (int i = 0; i < threadnum; i++) {
                    thpool_add_work(thpool, (void*) thread_func, (void*)&threadarr[i]);
                    
                }
                // let threads join (wait)
                thpool_wait(thpool);
                // free threads seq
                for (int i = 0; i < threadnum; i++) {
                    free(threadarr[i].obs_head);
                    free(threadarr[i].obs_seq);
                }
            }
            if (!feof(fp)) {
                threadarr[total % threadnum].obs_head = strdup(strtok(mystring, DELIMI));
                total++;
            }
            continue;
        } else {
            int bpcount = strlen(mystring);
            if(sequence_offset+bpcount>size) {
                size = sequence_offset+bpcount + 1;
                sequentie = realloc(sequentie,size);
            }
            while(mystring[bpcount - 1] == 10 || mystring[bpcount - 1]== 13) bpcount--;
            memcpy(sequentie + sequence_offset, mystring, bpcount);
            sequence_offset += bpcount;
            sequentie[sequence_offset] = '\0';
        }
    }
    thpool_destroy(thpool);
    fclose(fp);
    // print is uselless since all work is done, note total counts also the sequences that are of length 0!
    // printf("no. of seqs: %d\n", total);

    // close all file descriptors present in the data array of the threads. 
    for (int i = 0; i < threadnum; i++)  {
        fclose(threadarr[i].out);
        fclose(threadarr[i].aa);
        fclose(threadarr[i].dna);
    }
    if(threadnum > 1) combine(threadnum,out_header,threadarr);
    free(threadarr);
    printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (clock() - start) / (60.0 * CLOCKS_PER_SEC));
    return 0;
}

/**
* Function given to a thread during his creation. This thread will then execute this function.
*/
void* thread_func(void *threadarr) {
    thread_data *d = (thread_data*) threadarr;
    int cg = get_prob_from_cg(d->hmm, d->train, d->obs_seq);
    if (strlen(d->obs_seq)>70) viterbi(d->hmm, d->train, d->obs_seq, d->out, d->aa, d->dna, d->obs_head, d->wholegenome, cg, d->format);
    return NULL;
}

/**
* Error function:
* 1. Print error message
* 2. Call print_usage() from util_lib
* 3. EXIT program
*/
void print_error(const char* error_message, ...) {
    va_list list;
    va_start(list,error_message);
    vfprintf(stderr, error_message,list);
    va_end(list);
    print_usage();
    exit(EXIT_FAILURE);
}

/**
* Error function:
* 1. Print error message
* 2. EXIT program
*/
void print_file_error(const char* error_message, char* file) {
    fprintf(stderr, error_message,file);
    exit(EXIT_FAILURE);
}

/**
* This function will combine multiple files created by the invidual threads if number of threads > 1;
*/
void combine(int threadnum,char* out_header, thread_data *threadarr) {
    char aa_file[STRINGLEN] = "";
    char out_file[STRINGLEN] = "";
    char dna_file[STRINGLEN] = "";
    char mystring[STRINGLEN] = "";
    /* create output file name */
    strcpy(aa_file, out_header);
    strcat(aa_file, ".faa");
    strcpy(dna_file, out_header);
    strcat(dna_file, ".ffn");
    strcpy(out_file, out_header);
    strcat(out_file, ".out");

    remove (out_file);
    remove (aa_file);
    remove (dna_file);

    FILE * fp_aa = fopen (aa_file, "w");
    FILE * fp_out = fopen (out_file, "w");
    FILE * fp_dna = fopen (dna_file, "w");

    char ** lastline = malloc(threadnum* sizeof(char*));
    char ** currline =  malloc(threadnum* sizeof(char*));

    int j =0;

    for (int i = 0; i < threadnum; i++) {
        sprintf(mystring, "%s.out.tmp.%d", out_header, i);
        threadarr[i].out = fopen(mystring, "r");
        sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
        threadarr[i].aa = fopen(mystring, "r");
        sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
        threadarr[i].dna = fopen(mystring, "r");

        lastline[i] = malloc(STRINGLEN + 1);
        currline[i] = malloc(STRINGLEN + 1);
    }
    // Organize out file
    while (j != threadnum) {
        for (int i = 0; i < threadnum; i++) {
            if (lastline[i][0] != '\0') {
                fputs(lastline[i], fp_out);
                lastline[i][0] = '\0';
            }
            while(fgets(currline[i], STRINGLEN, threadarr[i].out)) {
                if (currline[i][0] == '>') {
                    memcpy(lastline[i], currline[i], strlen(currline[i]) + 1);
                    break;
                } else fputs(currline[i], fp_out);
            }
            if (feof(threadarr[i].out)) j++;
        }
    }
    // Organize faa file
    for (int i = 0; i < threadnum; i++) lastline[i][0] = '\0';
    j = 0;
    while (j != threadnum)    {
        for (int i = 0; i < threadnum; i++) {
            if (lastline[i][0] != '\0') {
                fputs(lastline[i], fp_aa);
                lastline[i][0] = '\0';
            }
            while(fgets(currline[i], STRINGLEN, threadarr[i].aa)) {
                if (currline[i][0] == '>') {
                    memcpy(lastline[i], currline[i], strlen(currline[i]) + 1);
                    break;
                }
                else {
                    fputs(currline[i], fp_aa);
                }
            }
            if (feof(threadarr[i].aa)) j++;
        }
    }

    // Organize dna file
    for (int i = 0; i < threadnum; i++) {
        lastline[i][0] = '\0';
    }
    j = 0;
    while (j != threadnum) {

        for (int i = 0; i < threadnum; i++) {
            if (lastline[i][0] != '\0') {
                fputs(lastline[i], fp_dna);
                lastline[i][0] = '\0';
            }
            while(fgets(currline[i], STRINGLEN, threadarr[i].dna)) {
                if (currline[i][0] == '>') {
                    memcpy(lastline[i], currline[i], strlen(currline[i]) + 1);
                    break;
                } else {
                    fputs(currline[i], fp_dna);
                }
            }
            if (feof(threadarr[i].dna))
            {
                j++;
            }
        }
    }

    //FREE
    for (int i = 0; i < threadnum; i++) {
        fclose(threadarr[i].out);
        fclose(threadarr[i].aa);
        fclose(threadarr[i].dna);
        sprintf(mystring, "%s.out.tmp.%d", out_header, i);
        sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
        sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
        free(lastline[i]);
        free(currline[i]);
    }

    fclose(fp_out);
    fclose(fp_aa);
    fclose(fp_dna);

    free(lastline);
    free(currline);


}