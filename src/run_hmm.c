#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "run_hmm.h"
#include "util_lib.h"

#include <pthread.h>

#define ADD_LEN 1024
#define STRINGLEN 4096+1
#define DELIMI " "
/**
* Entry point of program
* 1. Initialization of variables and datatypes
* 2. Check File acessiblity
* 3. call get_train_from_file, load hmm and train structswith data
* 4. read sequence file, start threads, run viterbi (in thread)
* 5. combine viterbi output files if nessecersairy.
*/
int main (int argc, char **argv) {
    /*Initialization datastructures*/
    TRAIN *train = malloc(sizeof(TRAIN));
    HMM *hmm = malloc(sizeof(HMM));

//    clock_t start = clock();

    /* number of threads */
    int threadnum = 1;

    char hmm_file[STRINGLEN*2] = "";
    char out_header[STRINGLEN-500] = "";

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
    char dna_file[STRINGLEN] = "";
    char out_file[STRINGLEN] = "";



    getcwd(train_dir,STRINGLEN);
    strcpy(mstate_file, train_dir);
    strcat(mstate_file, "/train/gene");
    strcpy(rstate_file, train_dir);
    strcat(rstate_file, "/train/rgene");
    strcpy(nstate_file, train_dir);
    strcat(nstate_file, "/train/noncoding");
    strcpy(sstate_file, train_dir);
    strcat(sstate_file, "/train/start");
    strcpy(pstate_file, train_dir);
    strcat(pstate_file, "/train/stop");
    strcpy(s1state_file, train_dir);
    strcat(s1state_file, "/train/stop1");
    strcpy(p1state_file, train_dir);
    strcat(p1state_file, "/train/start1");
    strcpy(dstate_file, train_dir);
    strcat(dstate_file, "/train/pwm");


    // arguments variables
    int wholegenome, format=0,d=0,e=0;

    // Get the command line variables
    int c;
    while ((c=getopt(argc, argv, "fd:e:w:t:p:")) != -1) {
        switch (c) {
        case 'w':
            wholegenome = atoi(optarg);
            if (wholegenome != 0 && wholegenome != 1)	print_error("ERROR: An incorrect value for the option -w was entered\n");
            break;
        case 'p':
            threadnum = atoi(optarg);
            if (threadnum < 1) print_error("ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
            //printf("Using %d threads.\n", threadnum);
            break;
        case 't':
            strcpy(train_file, optarg);
            snprintf(hmm_file,STRINGLEN*2,"%s/%s",train_dir,train_file);
            break;
        case 'f':
            format = 1;
            break;
        case 'e':
            e = 1;
            strcpy(out_file, optarg);

            break;
        case 'd':
            d = 1;
            strcpy(dna_file, optarg);
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

    // Allocate thread data structure
    thread_data *threadarr = malloc(sizeof(thread_data)* threadnum);

    // Initialize thread data structure
    for (int i = 0; i < threadnum; i++)   {
        if(threadnum > 1) {
            if(e) {
                sprintf(mystring, "%s.out.tmp.%d", out_header, i);
                threadarr[i].out = fopen(mystring, "w");
            }
            sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
            threadarr[i].aa = fopen(mystring, "w");
            if(d) {
                sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
                threadarr[i].dna = fopen(mystring, "w");
            }
        } else {
            if(e) {
                threadarr[i].out = fopen(out_file, "w");
            }
            sprintf(mystring, "%s.faa", out_header);
            threadarr[i].aa = stdout;
            if(d) {
                threadarr[i].dna = fopen(dna_file, "w");
            }
        }

        threadarr[i].hmm = hmm;
        threadarr[i].train = train;
        threadarr[i].d = d;
        threadarr[i].e = e;
        threadarr[i].wholegenome = wholegenome;
        threadarr[i].format = format;
    }


    // used for th sum of bpcounts
    int sequence_offset = 0;
    // totaal
    int  total=0;
    // array of threads (variables needed for pthreads)
    pthread_t thread[threadnum];
    // initial size of sequentie string, this will increase if there are strings in the file are bigger
    size_t size = 150;
    //sequentie
    char* sequentie = malloc(size* sizeof(char));

    while ( fgets (mystring, sizeof mystring, stdin)  ) {
        if (mystring[0] == '>'||feof(stdin)) {

            if (total > 0) {
                threadarr[(total-1) % threadnum].obs_seq = strdup(sequentie);
                sequence_offset = 0;
            }
            if ((total > 0 && total % threadnum == 0)||feof(stdin) ) {
                // Deal with the thread
                for (int i = 0; i < threadnum; i++) {
                    int rc = pthread_create(&thread[i], NULL, thread_func, (void*)&threadarr[i]);
                    if (rc) {
                        printf("Error: Unable to create thread, %d\n", rc);
                        exit(-1);
                    }
                }
                // let threads join (wait)
                for (int i = 0; i < threadnum; i++) {
                    int rc = pthread_join(thread[i], NULL);
                    if (rc) {
                        printf("Error: Unable to join threads, %d\n", rc);
                        exit(-1);
                    }
                }
                // free threads seq
                for (int i = 0; i < threadnum; i++) {
                    free(threadarr[i].obs_head);
                    free(threadarr[i].obs_seq);
                }
            }
            if(feof(stdin)) continue;
            threadarr[total % threadnum].obs_head = strdup(strtok(mystring, DELIMI));
            total++;
            continue;
        } else {
            //get length string
            int bpcount = strlen(mystring);
            // Als groter dan gealloceerd geheugen, maak het geheugen groter
            if(sequence_offset+bpcount>size) {
                size = sequence_offset+bpcount + 1;
                sequentie = realloc(sequentie,size);
            }
            //append
            while(mystring[bpcount - 1] == 10 || mystring[bpcount - 1]== 13) bpcount--;
            memcpy(sequentie + sequence_offset, mystring, bpcount);
            sequence_offset += bpcount;
            sequentie[sequence_offset] = '\0';
        }
    }
    free(sequentie);
    // print is uselless since all work is done, note total counts also the sequences that are of length 0!
    // printf("no. of seqs: %d\n", total);

    // close all file descriptors present in the data array of the threads.
    for (int i = 0; i < threadnum; i++)  {
        if (e) fclose(threadarr[i].out);
        fclose(threadarr[i].aa);
        if (d) fclose(threadarr[i].dna);
    }
    //combine output files if nesscairy
    if(threadnum > 1) combine(threadnum,out_header,dna_file,out_file,threadarr,d,e);
    free(threadarr);
    free(hmm);
    free(train);
    //printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (clock() - start) / (60.0 * CLOCKS_PER_SEC));
    return 0;
}

/**
* Function given to a thread during his creation. This thread will then execute this function.
*/
void* thread_func(void *threadarr) {
    thread_data *d = (thread_data*) threadarr;
    int cg = get_prob_from_cg(d->hmm, d->train, d->obs_seq);
    if (strlen(d->obs_seq)>70) viterbi(d->hmm, d->train, d->obs_seq, d->out, d->aa, d->dna, d->obs_head, d->wholegenome, cg, d->format,d->d,d->e);
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
void combine(int threadnum,char* out_header,char* dna_file,char* out_file, thread_data *threadarr,int d,int e) {
    char aa_file[STRINGLEN] = "";
    char mystring[STRINGLEN] = "";
    /* create output file name */
    strcpy(aa_file, out_header);
    strcat(aa_file, ".faa");

    FILE *fp_dna,*fp_out;
    if (d) {
        fp_dna = fopen (dna_file, "w");
    }
    if(e) {
        fp_out = fopen (out_file, "w");
    }

    char ** lastline = malloc(threadnum* sizeof(char*));
    char * currline =  calloc(STRINGLEN, sizeof(char));

    int j =0;

    for (int i = 0; i < threadnum; i++) {
        if(e) {
            sprintf(mystring, "%s.out.tmp.%d", out_header, i);
            threadarr[i].out = fopen(mystring, "r");
        }
        sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
        threadarr[i].aa = fopen(mystring, "r");
        if(d) {
            sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
            threadarr[i].dna = fopen(mystring, "r");
        }
        lastline[i] = calloc(STRINGLEN, sizeof(char));
    }
    FILE * fp_aa = fopen (aa_file, "w");
    // Organize out file
    if(e) {
        while (j != threadnum) {
            for (int i = 0; i < threadnum; i++) {
                if (lastline[i][0] != '\0') {
                    fputs(lastline[i], fp_out);
                    lastline[i][0] = '\0';
                }
                while(fgets(currline, STRINGLEN, threadarr[i].out)) {
                    if (currline[0] == '>') {
                        memcpy(lastline[i], currline, strlen(currline) + 1);
                        break;
                    } else {
                        fputs(currline, fp_out);
                    }
                }
                if (feof(threadarr[i].out)) j++;
            }
        }
    }
    // Organize faa file
    for (int i = 0; i < threadnum; i++) lastline[i][0] = '\0';
    j = 0;
    while (j != threadnum)    {
        for (int i = 0; i < threadnum; i++) {
            if (lastline[i][0] != '\0') {
                printf("%s\n",lastline[i]);
                lastline[i][0] = '\0';
            }
            while(fgets(currline, STRINGLEN, threadarr[i].aa)) {
                if (currline[0] == '>') {
                    memcpy(lastline[i], currline, strlen(currline) + 1);
                    break;
                }
                else {
                    printf("%s\n",currline);
                }
            }
            if (feof(threadarr[i].aa)) j++;
        }
    }

    // Organize dna file
    if(d) {
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
                while(fgets(currline, STRINGLEN, threadarr[i].dna)) {
                    if (currline[0] == '>') {
                        memcpy(lastline[i], currline, strlen(currline) + 1);
                        break;
                    } else {
                        fputs(currline, fp_dna);
                    }
                }
                if (feof(threadarr[i].dna)) j++;
            }
        }
    }
    //FREE
    for (int i = 0; i < threadnum; i++) {
        if (e) {
            fclose(threadarr[i].out);
            fclose(threadarr[i].aa);
            sprintf(mystring, "%s.out.tmp.%d", out_header, i);
            remove(mystring);
        }
        if (d) {
            fclose(threadarr[i].dna);
            sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
            remove(mystring);
        }
        sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
        remove(mystring);
        free(lastline[i]);
    }

    if(e)    fclose(fp_out);
    fclose(fp_aa);
    if(d)    fclose(fp_dna);
    free(lastline);
    free(currline);
}

