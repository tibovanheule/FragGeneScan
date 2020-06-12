#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "run_hmm.h"
#include "util_lib.h"

#include <pthread.h>

#define ADD_LEN 1024
#define STRINGLEN 4096

/**
* Entry point of program
* 1. Initialization of variables and datatypes
* 2. Check File acessiblity
*/
int main (int argc, char **argv) {
    //file handlers
    FILE *fp_out, *fp_aa, *fp_dna, *fp;
    /*Initialization datastructures*/
    TRAIN train;
    HMM hmm;
    thread_data *threadarr;

    clock_t start = clock();
    int i, j, max;
    int wholegenome;
    int format=0;
    int count=0;
    int currcount = 0;
    int total = 0;
    int *obs_seq_len;
    /* count the length of each line in input file */
    int bp_count;
    /* number of threads */
    int threadnum = 1;
    char *obs_seq, *obs_head;
    char hmm_file[STRINGLEN] = "";
    char out_header[STRINGLEN] = "";
    char aa_file[STRINGLEN] = "";
    char seq_file[STRINGLEN] = "";
    char out_file[STRINGLEN] = "";
    char dna_file[STRINGLEN] = "";
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
    char **lastline, **currline;

    strncpy(train_dir, argv[0], strlen(argv[0])-12);
    strcat(train_dir, "train/");
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
            if (access(hmm_file, F_OK)==-1) print_error("ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
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
    hmm.N=NUM_STATE;
    /* Read the given files and give HMM and TRAIN */
    get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

    // Initialize thread data structure
    threadarr = (thread_data*) calloc(sizeof(thread_data), threadnum);
    memset(threadarr, '\0', sizeof(thread_data) * threadnum);
    for (i = 0; i < threadnum; i++)   {
        if(threadnum > 1) sprintf(mystring, "%s.out.tmp.%d", out_header, i);
        else sprintf(mystring, "%s.out", out_header);

        threadarr[i].out = fopen(mystring, "w");

        if(threadnum > 1) sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
        else sprintf(mystring, "%s.faa", out_header);

        threadarr[i].aa = fopen(mystring, "w");

        if(threadnum > 1) sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
        else sprintf(mystring, "%s.ffn", out_header);

        threadarr[i].dna = fopen(mystring, "w");

        threadarr[i].hmm = (HMM*) malloc(sizeof(HMM));
        threadarr[i].hmm = &hmm;
        threadarr[i].train = (TRAIN*) malloc(sizeof(TRAIN));
        threadarr[i].train = &train;

        threadarr[i].wholegenome = wholegenome;
        threadarr[i].format = format;
    }

    pthread_t thread[threadnum];

    fp = fopen (seq_file, "r");
    // tel gewoon aantal reads
    while ( fgets (mystring, sizeof mystring, fp) ) {
        if (mystring[0] == '>')      count++;
    }


    // lengte van elke sequentie
    obs_seq_len = (int *) malloc(count * sizeof(int));
    printf("no. of seqs: %d\n", count);
    // TODO controle

    i = 0;
    count = 0;
    rewind(fp);
    while ( fgets (mystring, sizeof mystring, fp) ) {
        if (mystring[0] == '>') {
            if (i>0) {
                obs_seq_len[count] = i;
                count++;
                i = 0;
            }

        } else {
            bp_count = strlen(mystring);
            while(mystring[bp_count-1] == 10 || mystring[bp_count-1]==13)	bp_count --;
            i += bp_count;
        }
    }
    obs_seq_len[count] = i;

    rewind(fp);
    total = 0;
    count = 0;
    j = 0;

    while (!(feof(fp))) {
        memset(mystring, '\0', sizeof mystring);
        fgets (mystring, sizeof mystring, fp);
        bp_count = strlen(mystring);
        while(mystring[bp_count - 1] == 10 || mystring[bp_count - 1]== 13) bp_count--;

        if (mystring[0] == '>' || feof(fp)) {
            if (feof(fp)) {
                // append to seq
                memcpy(threadarr[currcount].obs_seq + j, mystring, bp_count);
                j += bp_count;
            }
            if ((count > 0 && count % threadnum == 0) || feof(fp)) {
                // Deal with the thread
                for (i = 0; i < count; i++) {
                    int rc = pthread_create(&thread[i], NULL, thread_func, (void*)&threadarr[i]);
                    // TODO print error with custom error code
                    if (rc) {
                        printf("Error: Unable to create thread, %d\n", rc);
                        exit(-1);
                    }
                }
                // let threads join (wait)
                for (i = 0; i < count; i++) {
                    int rc = pthread_join(thread[i], NULL);
                    if (rc) {
                        printf("Error: Unable to join threads, %d\n", rc);
                        exit(-1);
                    }
                }
                // free threads array
                for (i = 0; i < count; i++) {
                    free(threadarr[i].obs_head);
                    free(threadarr[i].obs_seq);
                    threadarr[i].obs_head = NULL;
                    threadarr[i].obs_seq = NULL;
                }
                count = 0;
            }
            if (!(feof(fp))) {
                threadarr[count].obs_head = (char *) malloc((bp_count+1)* sizeof(char));
                memcpy(threadarr[count].obs_head, mystring, bp_count);
                //threadarr[count].obs_seq = NULL;
                threadarr[count].obs_seq = (char*) calloc((obs_seq_len[total] + 1), sizeof(char));
                total++;
                currcount = count;
                count++;
                j = 0;
                max = 0;
            }

        } else {
            // append aan de sequentie
            memcpy(threadarr[currcount].obs_seq + j, mystring, bp_count);
            j += bp_count;
        }
        if(feof(fp)) break;
    }
    //if (threadnum == 1) fclose(fp);

    for (i = 0; i < threadnum; i++)  {
        fclose(threadarr[i].out);
        fclose(threadarr[i].aa);
        fclose(threadarr[i].dna);
    }


    //uitschrijven
    if(threadnum > 1) {
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

        fp_aa = fopen (aa_file, "w");
        fp_out = fopen (out_file, "w");
        fp_dna = fopen (dna_file, "w");

        lastline = (char**) malloc(threadnum* sizeof(char*));
        memset(lastline, '\0', sizeof(char*) * threadnum);
        currline = (char**)malloc(threadnum* sizeof(char*));
        memset(currline, '\0', sizeof(char*) * threadnum);
        for (i = 0; i < threadnum; i++) {
            sprintf(mystring, "%s.out.tmp.%d", out_header, i);
            threadarr[i].out = fopen(mystring, "r");
            sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
            threadarr[i].aa = fopen(mystring, "r");
            sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
            threadarr[i].dna = fopen(mystring, "r");

            lastline[i] = (char*) calloc(STRINGLEN + 1, sizeof(char));
            currline[i] = (char*) calloc(STRINGLEN + 1, sizeof(char));
        }

        // Organize out file
        while (1) {
            j = 0;
            for (i = 0; i < threadnum; i++) {
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
            if (j == threadnum) {
                break;
            }
        }
        // Organize faa file
        for (i = 0; i < threadnum; i++) lastline[i][0] = '\0';
        while (1)    {
            j = 0;
            for (i = 0; i < threadnum; i++) {
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
            if (j == threadnum)
            {
                break;
            }
        }

        // Organize dna file
        for (i = 0; i < threadnum; i++) {
            lastline[i][0] = '\0';
        }
        while (j != threadnum) {
            j = 0;
            for (i = 0; i < threadnum; i++) {
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

        for (i = 0; i < threadnum; i++) {
            fclose(threadarr[i].out);
            fclose(threadarr[i].aa);
            fclose(threadarr[i].dna);
            sprintf(mystring, "%s.out.tmp.%d", out_header, i);
            remove(mystring);
            sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
            remove(mystring);
            sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
            remove(mystring);
            free(lastline[i]);
            free(currline[i]);
        }
        free(threadarr);
        free(lastline);
        free(currline);

        free(obs_seq_len);
        //free(obs_head);
        fclose(fp_out);
        fclose(fp_aa);
        fclose(fp_dna);
    }
    fclose(fp);
    printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (clock() - start) / (60.0 * CLOCKS_PER_SEC));
    return 0;
}


void* thread_func(void *threadarr) {
    thread_data *d = (thread_data*) threadarr;
    d->cg = get_prob_from_cg(d->hmm, d->train, d->obs_seq); /* cg - 26 Ye April 16, 2016 */
    if (strlen(d->obs_seq)>70) viterbi(d->hmm, d->train, d->obs_seq, d->out, d->aa, d->dna, d->obs_head, d->wholegenome, d->cg, d->format);
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