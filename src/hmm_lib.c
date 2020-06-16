#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "hmm.h"
#include "util_lib.h"
#include "hmm_lib.h"

const double log53 = -0.634878; //log(0.53);
const double log16 = -1.832581; //log(0.16);
const double log30 = -1.203973; //log(0.30);
const double log25 = -1.386294; //log(0.25);
const double log95 = -0.051293; //log(0.95);
const double log54 = -0.616186; //log(0.54);
const double log83 = -0.186330; //log(0.83);
const double log07 = -2.659260; //log(0.07);

void viterbi(HMM *hmm_ptr, TRAIN *train_ptr, char *O, FILE *fp_out, FILE *fp_aa, FILE *fp_dna, char *head, int whole_genome, int cg, int format) {
    double max_dbl = INFINITY;
    int len_seq = strlen(O);
    double ** alpha = dmatrix(NUM_STATE, len_seq);                      /* viterbi prob array */
    int ** path = imatrix(NUM_STATE, len_seq);                          /* viterbi path array */
    int j;
    int temp_t;
    double prob_save, start_freq;
    int from, from0, to;   /*from0: i-2 position, from: i-1 position */
    int from2;             /* from2: i-2, i-1 for condition in emission probability */
    double temp_alpha;
    int gene_len = 60;
    int num_d;          /* the number of delete */
    int freq_id;
    double h_kd, r_kd, p_kd;


    int start_t, dna_start_t, dna_start_t_withstop;
    int end_t, dna_end_t;
    int prev_match;
    int start_orf;
    double final_score;

    int insert[100];
    int delete[100];
    int insert_id, delete_id;

    int temp_i[6] = {0,0,0,0,0,0};
    int temp_i_1[6] = {0,0,0,0,0,0};

    int num_N=0;
    /***************************************************************/
    /* initialize                                                  */
    /***************************************************************/
    int refine = 0;

    if (whole_genome==1) {
        gene_len = 120;
        refine = 1;
    }

    for (int i=0; i<NUM_STATE; i++)         alpha[i][0] = -1 * hmm_ptr->pi[i];

    /* stop state */
    if ((O[0] == 'T'|| O[0] == 't')  &&
            (((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'A'|| O[2] == 'a')) ||
             ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'G'|| O[2] == 'g')) ||
             ((O[1] == 'G'|| O[1] == 'g') && (O[2] == 'A'|| O[2] == 'a')))) {

        alpha[E_STATE][0] = INFINITY;
        alpha[E_STATE][1] = INFINITY;
        path[E_STATE][1] = E_STATE;
        path[E_STATE][2] = E_STATE;

        alpha[M6_STATE][2] = INFINITY;
        alpha[M5_STATE][1] = INFINITY;
        alpha[M4_STATE][0] = INFINITY;
        alpha[M3_STATE][2] = INFINITY;
        alpha[M2_STATE][1] = INFINITY;
        alpha[M1_STATE][0] = INFINITY;

        if ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'A'|| O[2] == 'a')) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - log53;
        } else if ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'G'|| O[2] == 'g')) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - log16;
        } else if((O[1] == 'G'|| O[1] == 'g') && (O[2] == 'A'|| O[2] == 'a')) {
            alpha[E_STATE][2] = alpha[E_STATE][2] - log30;
        }
    }

    if ((O[2] == 'A'|| O[0] == 'a')  &&
            (((O[0] == 'T'|| O[0] == 't') && (O[1] == 'T'|| O[1] == 't')) ||
             ((O[0] == 'C'|| O[0] == 'c') && (O[1] == 'T'|| O[1] == 't')) ||
             ((O[0] == 'T'|| O[0] == 't') && (O[1] == 'C'|| O[1] == 'c')))) {

        alpha[S_STATE_1][0] = INFINITY;
        alpha[S_STATE_1][1] = INFINITY;
        alpha[S_STATE_1][2] = alpha[S_STATE][0];
        path[S_STATE_1][1] = S_STATE_1;
        path[S_STATE_1][2] = S_STATE_1;

        alpha[M3_STATE_1][2] = INFINITY;
        alpha[M6_STATE_1][2] = INFINITY;

        if ((O[0] == 'T'|| O[0] == 't') && (O[1] == 'T'|| O[1] == 't')) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log53;
        } else if ((O[0] == 'C'|| O[0] == 'c') && (O[1] == 'T'|| O[1] == 't')) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log16;
        } else if((O[0] == 'T'|| O[0] == 't') && (O[1] == 'C'|| O[1] == 'c')) {
            alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log30;
        }
    }

    /******************************************************************/
    /*  fill out the rest of the columns                              */
    /******************************************************************/
    int *vpath = ivector(len_seq);/* optimal path after backtracking */
    for (int t = 1; t < len_seq; t++) {
        from = nt2int(O[t-1]);
        from0 = (t>1)? nt2int(O[t-2]):2;
        to = nt2int(O[t]);

        /* if DNA is other than ACGT, do it later */
        if (from==4)          from=2;

        if (from0==4)             from0=2;

        if (to==4) {
            to=2;
            num_N += 1;
        } else {
            num_N=0;
        }
        from2 = from0*4+from;

        /******************/
        /* M state        */
        /******************/

        for (int i=M1_STATE; i<=M6_STATE; i++)   {

            if (alpha[i][t]<max_dbl) {

                if (t!=0) {

                    if (i==M1_STATE) {

                        /* from M state */
                        j = M6_STATE;
                        alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - train_ptr->trans[cg][0][from2][to];
                        path[i][t] = j;

                        /* from D state */
                        if (whole_genome==0) {
                            for (int j=M5_STATE; j>=M1_STATE; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1<i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if(num_d>0) {
                                    temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - train_ptr->trans[cg][0][from2][to]
                                                 - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }

                        /* from Start state */
                        temp_alpha = alpha[S_STATE][t-1] - train_ptr->trans[cg][0][from2][to];
                        if ( temp_alpha < alpha[i][t] ) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = S_STATE;
                        }

                    } else {  /*i ==M2-M6*/

                        /* from M state */
                        j = i - 1;
                        alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_MM] - train_ptr->trans[cg][i-M1_STATE][from2][to];
                        path[i][t] = j;


                        /* from D state */
                        if (whole_genome==0) {
                            for (int j=M6_STATE; j>=M1_STATE; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1 < i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d > 0) {
                                    temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - train_ptr->trans[cg][i-M1_STATE][from2][to]
                                                 - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }
                    }

                    /* from I state */
                    j = (i==M1_STATE)? I6_STATE:I1_STATE + (i - M1_STATE -1);

                    /* to aviod stop codon */
                    if (t<2) {
                    } else if((i==M2_STATE || i==M5_STATE) && (O[temp_i[j-I1_STATE]] == 'T'||O[temp_i[j-I1_STATE]] =='t') &&
                              (((O[t] == 'A'||O[t] == 'a') && (O[t+1] =='A'||O[t+1] =='a')) ||
                               ((O[t] == 'A'||O[t] == 'a') && (O[t+1] =='G'||O[t+1] =='g')) ||
                               ((O[t] == 'G'||O[t] == 'g') && (O[t+1] =='A'||O[t+1] =='a')))) {

                    } else if ((i==M3_STATE || i==M6_STATE) && (O[temp_i[j-I1_STATE]-1] == 'T'||O[temp_i[j-I1_STATE]-1] =='t') &&
                               (((O[temp_i[j-I1_STATE]] == 'A'||O[temp_i[j-I1_STATE]] == 'a') && (O[t] =='A'||O[t] == 'a')) ||
                                ((O[temp_i[j-I1_STATE]] == 'A'||O[temp_i[j-I1_STATE]] == 'a') && (O[t] =='G'||O[t] == 'g')) ||
                                ((O[temp_i[j-I1_STATE]] == 'G'||O[temp_i[j-I1_STATE]] == 'g') && (O[t] =='A'||O[t] == 'a')))) {
                    } else {
                        temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_IM] - log25;
                        if ( temp_alpha < alpha[i][t]) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;
                        }
                    }
                }
            }
        }

        /******************/
        /* I state        */
        /******************/
        for (int i=I1_STATE; i<=I6_STATE; i++) {

            if (t!=0) {
                /* from I state */
                j = i;
                alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
                path[i][t] = j;

                /* from M state */
                j = i - I1_STATE + M1_STATE ;
                if (i==I6_STATE) {
                    temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                } else {
                    temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
                }
                if (temp_alpha < alpha[i][t]) {
                    alpha[i][t] = temp_alpha;
                    path[i][t] = j;

                    temp_i[i-I1_STATE] = t-1;
                }
            }
        }

        /******************/
        /* M' state        */
        /******************/

        for (int i=M1_STATE_1; i<=M6_STATE_1; i++)   {
            if  ((i==M1_STATE_1 || i==M4_STATE_1)&& t>=3 &&
                    (((O[t-3] == 'T'||O[t-3] == 't') && (O[t-2] == 'T'||O[t-2] == 't') && (O[t-1] == 'A'||O[t-1] =='a')) ||
                     ((O[t-3] == 'C'||O[t-3] == 'c') && (O[t-2] == 'T'||O[t-2] == 't') && (O[t-1] == 'A'||O[t-1] =='a')) ||
                     ((O[t-3] == 'T'||O[t-3] == 't') && (O[t-2] == 'C'||O[t-2] == 'c') && (O[t-1] == 'A'||O[t-1] =='a')))) {

                /* from Start state  since this is actually stop codon in minus strand */
                alpha[i][t] = alpha[S_STATE_1][t-1] - train_ptr->rtrans[cg][i-M1_STATE_1][from2][to];
                path[i][t] = S_STATE_1;

            } else {

                if (t!=0) {
                    if (i==M1_STATE_1 ) {

                        /* from M state */
                        j = M6_STATE_1;
                        alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - train_ptr->rtrans[cg][0][from2][to];
                        path[i][t] = j;

                        /* from D state */
                        if (whole_genome==0) {
                            for (int j=M5_STATE_1; j>=M1_STATE_1; j--) {
                                if (j >= i) {
                                    num_d = i-j+6;
                                } else if (j+1 <i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d > 0) {
                                    temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - train_ptr->rtrans[cg][0][from2][to]
                                                 - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }

                    } else {

                        /* from M state */
                        j = i - 1;
                        alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_MM] - train_ptr->rtrans[cg][i-M1_STATE_1][from2][to];
                        path[i][t] = j;

                        /* from D state */
                        if (whole_genome==0) {
                            for (int j=M6_STATE_1; j>=M1_STATE_1; j--) {
                                if (j >= i ) {
                                    num_d = i-j+6;
                                } else if (j+1 < i) {
                                    num_d = i-j;
                                } else {
                                    num_d = -10;
                                }
                                if (num_d>0) {
                                    temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - train_ptr->rtrans[cg][i-M1_STATE_1][from2][to]
                                                 - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
                                    if ( temp_alpha < alpha[i][t]) {
                                        alpha[i][t] = temp_alpha;
                                        path[i][t] = j;
                                    }
                                }
                            }
                        }
                    }

                    /* from I state */
                    if (i==M1_STATE_1) {
                        j = I6_STATE_1;
                    } else {
                        j = I1_STATE_1 + (i - M1_STATE_1 -1);
                    }


                    /* to aviod stop codon */
                    /* to aviod stop codon */
                    if (t<2) {
                    } else  if((i==M2_STATE_1 || i==M5_STATE_1) && (O[t+1] == 'A'||O[t+1] == 'a') &&
                               (((O[temp_i_1[j-I1_STATE_1]] == 'T'|| O[temp_i_1[j-I1_STATE_1]] == 't') && (O[t] =='T'|| O[t] =='t')) ||
                                ((O[temp_i_1[j-I1_STATE_1]] == 'C'|| O[temp_i_1[j-I1_STATE_1]] == 'c') && (O[t] =='T'|| O[t] =='t')) ||
                                ((O[temp_i_1[j-I1_STATE_1]] == 'T'|| O[temp_i_1[j-I1_STATE_1]] == 't') && (O[t] =='C'|| O[t] =='c')))) {

                    } else if ((i==M3_STATE_1 || i==M6_STATE_1) && (O[t] == 'A'||O[t] == 'a') &&
                               (((O[temp_i_1[j-I1_STATE_1]-1] == 'T'|| O[temp_i_1[j-I1_STATE_1]-1]=='t') &&
                                 (O[temp_i_1[j-I1_STATE_1]] =='T'|| O[temp_i_1[j-I1_STATE_1]] =='t')) ||
                                ((O[temp_i_1[j-I1_STATE_1]-1] == 'C'|| O[temp_i_1[j-I1_STATE_1]-1]=='c') &&
                                 (O[temp_i_1[j-I1_STATE_1]] =='T'|| O[temp_i_1[j-I1_STATE_1]] =='t')) ||
                                ((O[temp_i_1[j-I1_STATE_1]-1] == 'T'|| O[temp_i_1[j-I1_STATE_1]-1]=='t') &&
                                 (O[temp_i_1[j-I1_STATE_1]] =='C'|| O[temp_i_1[j-I1_STATE_1]] =='c')))) {
                    } else {

                        temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_IM] - log25;
                        if ( temp_alpha < alpha[i][t]) {
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;
                        }
                    }
                }
            }
        }

        /******************/
        /* I' state        */
        /******************/
        for (int i=I1_STATE_1; i<=I6_STATE_1; i++) if (t!=0) {
                /* from I state */
                j = i;
                alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
                path[i][t] = j;

                /* from M state */
                if (path[S_STATE_1][t-3]!= R_STATE && path[S_STATE_1][t-4] !=R_STATE && path[S_STATE_1][t-5] !=R_STATE) {
                    j = i - I1_STATE_1 + M1_STATE_1;
                    if (i==I6_STATE_1) {
                        temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                    } else {
                        temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                    }
                    if (temp_alpha < alpha[i][t]) {
                        alpha[i][t] = temp_alpha;
                        path[i][t] = j;

                        temp_i_1[i-I1_STATE_1] = t-1;
                    }
                }

            }

        /***********************/
        /* Non_coding state    */
        /***********************/

        if (t!=0) {
            alpha[R_STATE][t] = alpha[R_STATE][t-1] - train_ptr->noncoding[cg][from][to] - hmm_ptr->tr[TR_RR];
            path[R_STATE][t] = R_STATE;

            temp_alpha = alpha[E_STATE][t-1]  - hmm_ptr->tr[TR_ER];
            if (temp_alpha < alpha[R_STATE][t] ) {
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE;
            }

            temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ER] ;
            if (temp_alpha < alpha[R_STATE][t] ) {
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE_1;
            }
            alpha[R_STATE][t] -= log95;
        }

        /******************/
        /* END state      */
        /******************/
        if (alpha[E_STATE][t] == 0) {

            alpha[E_STATE][t] = max_dbl;
            path[E_STATE][t] = NOSTATE;

            if (t < len_seq -2 && (O[t] == 'T'||O[t] == 't')  &&
                    (((O[t+1] == 'A'||O[t+1] == 'a') && (O[t+2] == 'A'||O[t+2] =='a')) ||
                     ((O[t+1] == 'A'||O[t+1] == 'a') && (O[t+2] == 'G'||O[t+2] =='g')) ||
                     ((O[t+1] == 'G'||O[t+1] == 'g') && (O[t+2] == 'A'||O[t+2] =='a')))) {

                alpha[E_STATE][t+2] = max_dbl;
                /* transition from frame4,frame5,and frame6 */
                temp_alpha = alpha[M6_STATE][t-1] - hmm_ptr->tr[TR_GE];
                if (temp_alpha < alpha[E_STATE][t+2]) {
                    alpha[E_STATE][t+2] = temp_alpha;
                    path[E_STATE][t] = M6_STATE;
                }

                /* transition from frame1,frame2,and frame3 */
                temp_alpha  = alpha[M3_STATE][t-1] - hmm_ptr->tr[TR_GE];
                if (temp_alpha < alpha[E_STATE][t+2]) {
                    alpha[E_STATE][t+2] = temp_alpha;
                    path[E_STATE][t] = M3_STATE;
                }

                alpha[E_STATE][t] = max_dbl;
                alpha[E_STATE][t+1] = max_dbl;
                path[E_STATE][t+1] = E_STATE;
                path[E_STATE][t+2] = E_STATE;

                alpha[M6_STATE][t+2] = max_dbl;
                alpha[M5_STATE][t+1] = max_dbl;
                alpha[M4_STATE][t] = max_dbl;
                alpha[M3_STATE][t+2] = max_dbl;
                alpha[M2_STATE][t+1] = max_dbl;
                alpha[M1_STATE][t] = max_dbl;

                if ((O[t+1] == 'A'||O[t+1] =='a') && (O[t+2] == 'A'||O[t+2] =='a')) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log54;
                } else if ((O[t+1] == 'A'||O[t+1] =='a') && (O[t+2] == 'G'||O[t+2] =='g')) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log16;
                } else if((O[t+1] == 'G'||O[t+1] == 'g') && (O[t+2] == 'A'||O[t+2] =='a')) {
                    alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log30;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;

                double sub_sum = 0;
                int sub_count = 0;

                if (t>=60) { /* bug reported by Yu-Wei */
                    for(int i=-60; i<=-3; i++) {
                        if (t+i+2 < len_seq)
                        {
                            start_freq -= train_ptr->stop[cg][i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                        }
                    }
                } else {
                    for(int i=(-1*t); i<=-3; i++) {
                        if (t+i+2 < len_seq)
                        {
                            sub_sum += train_ptr->stop[cg][i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                        }
                    }
                    sub_sum = sub_sum * 58 / (-3 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = train_ptr->E_dist[cg][2] * exp(-1*pow(start_freq-train_ptr->E_dist[cg][1],2)/(2*pow(train_ptr->E_dist[cg][0],2)));
                r_kd = train_ptr->E_dist[cg][5] * exp(-1*pow(start_freq-train_ptr->E_dist[cg][4],2)/(2*pow(train_ptr->E_dist[cg][3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log(p_kd);
            }
        }

        /*************************************************/
        /* START' state                                  */
        /* origianlly stop codon of genes in - strand    */
        /*************************************************/
        if (alpha[S_STATE_1][t] == 0) {

            alpha[S_STATE_1][t] = max_dbl;
            path[S_STATE_1][t] = NOSTATE;


            if (t<len_seq-2 && (O[t+2] == 'A'||O[t+2] == 'a') &&
                    (((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'T'|| O[t+1] == 't')) ||
                     ((O[t] == 'C'||O[t] =='c') && (O[t+1] == 'T'|| O[t+1] == 't')) ||
                     ((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'C'|| O[t+1] == 'c')))) {

                alpha[S_STATE_1][t] = max_dbl;
                path[S_STATE_1][t] = R_STATE;
                alpha[S_STATE_1][t+1] = max_dbl;
                alpha[S_STATE_1][t+2] = alpha[R_STATE][t-1] - hmm_ptr->tr[TR_RS];
                path[S_STATE_1][t+1] = S_STATE_1;
                path[S_STATE_1][t+2] = S_STATE_1;

                temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ES];
                if (temp_alpha < alpha[S_STATE_1][t+2]) {
                    alpha[S_STATE_1][t+2] = temp_alpha;
                    path[S_STATE_1][t] = E_STATE_1;
                }

                temp_alpha = alpha[E_STATE][t-1] - hmm_ptr->tr[TR_ES1];
                if (temp_alpha < alpha[S_STATE_1][t+2]) {
                    alpha[S_STATE_1][t+2] = temp_alpha;
                    path[S_STATE_1][t] = E_STATE;
                }

                alpha[M3_STATE_1][t+2] = max_dbl;
                alpha[M6_STATE_1][t+2] = max_dbl;

                if ((O[t] == 'T'||O[t] == 't') && (O[t+1] == 'T'||O[t+1] == 't')) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log54;
                } else if ((O[t] == 'C'||O[t] =='c') && (O[t+1] == 'T'||O[t+1]=='t')) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log16;
                } else if((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'C'||O[t+1] =='c')) {
                    alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log30;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;
                for(int i=3; i<=60; i++) {
                    if (t+i+2 < len_seq)
                    {
                        start_freq -= train_ptr->start1[cg][i-3][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                    }
                }
                h_kd = train_ptr->S1_dist[cg][2] * exp(-1*pow(start_freq-train_ptr->S1_dist[cg][1],2)/(2*pow(train_ptr->S1_dist[cg][0],2)));
                r_kd = train_ptr->S1_dist[cg][5] * exp(-1*pow(start_freq-train_ptr->S1_dist[cg][4],2)/(2*pow(train_ptr->S1_dist[cg][3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log(p_kd);
            }
        }

        /************************/
        /* START state          */
        /************************/
        if (alpha[S_STATE][t] == 0) {

            alpha[S_STATE][t] = max_dbl;
            path[S_STATE][t] = NOSTATE;

            if (t<len_seq-2 &&  (O[t+1] == 'T'||O[t+1] =='t') && (O[t+2] == 'G'||O[t+2] =='g')&&
                    ((O[t] == 'A'||O[t] =='a') || (O[t] == 'G'||O[t] =='g') ||  (O[t] == 'T'||O[t] =='t'))) {

                alpha[S_STATE][t] = max_dbl;
                alpha[S_STATE][t+1] = max_dbl;
                alpha[S_STATE][t+2] = alpha[R_STATE][t-1] - hmm_ptr->tr[TR_RS];
                path[S_STATE][t] = R_STATE;
                path[S_STATE][t+1] = S_STATE;
                path[S_STATE][t+2] = S_STATE;

                temp_alpha = alpha[E_STATE][t-1] - hmm_ptr->tr[TR_ES];
                if (temp_alpha < alpha[S_STATE][t+2]) {
                    alpha[S_STATE][t+2] = temp_alpha;
                    path[S_STATE][t] = E_STATE;
                }

                temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ES1];
                if (temp_alpha < alpha[S_STATE][t+2]) {
                    alpha[S_STATE][t+2] = temp_alpha;
                    path[S_STATE][t] = E_STATE_1;
                }


                if ((O[t] == 'A'||O[t] =='a') ) {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log83;
                } else if ((O[t] == 'G'||O[t] =='g')) {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(0.10);
                } else if((O[t] == 'T'||O[t] == 't')) {
                    alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log07;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;
                double sub_sum = 0;
                int sub_count = 0;

                if (t>=30) {
                    for(int i=-30; i<=30; i++) {
                        if (t+i+2 < len_seq)
                        {
                            start_freq -= train_ptr->start[cg][i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                        }
                    }
                } else {
                    for(int i=(-1*t); i<=30; i++) {
                        if (t+i+2 < len_seq)
                        {
                            sub_sum += train_ptr->start[cg][i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                        }
                    }
                    sub_sum = sub_sum * 61 / (30 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = train_ptr->S_dist[cg][2] * exp(-1*pow(start_freq-train_ptr->S_dist[cg][1],2)/(2*pow(train_ptr->S_dist[cg][0],2)));
                r_kd = train_ptr->S_dist[cg][5] * exp(-1*pow(start_freq-train_ptr->S_dist[cg][4],2)/(2*pow(train_ptr->S_dist[cg][3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(p_kd);

            }
        }

        /**********************************************/
        /* END' state                                 */
        /* originally start codon of genes in - strand */
        /**********************************************/
        if (alpha[E_STATE_1][t] == 0) {

            alpha[E_STATE_1][t] = max_dbl;
            path[E_STATE_1][t] = NOSTATE;

            if (t < len_seq - 2 && (O[t] == 'C'||O[t] =='c') && (O[t+1] == 'A'||O[t+1] == 'a') &&
                    ((O[t+2] == 'T'||O[t+2] =='t') || (O[t+2] == 'C'||O[t+2] =='c') || (O[t+2] == 'A'||O[t+2] =='a'))) {

                /* transition from frame6 */
                alpha[E_STATE_1][t+2] = alpha[M6_STATE_1][t-1] - hmm_ptr->tr[TR_GE];
                path[E_STATE_1][t] = M6_STATE_1;
                alpha[E_STATE_1][t] = max_dbl;
                alpha[E_STATE_1][t+1] = max_dbl;
                path[E_STATE_1][t+1] = E_STATE_1;
                path[E_STATE_1][t+2] = E_STATE_1;

                if ((O[t+2] == 'T'||O[t+2] == 't') ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log83;
                } else if ((O[t+2] == 'C'||O[t+2] =='c') ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(0.10);
                } else if((O[t+2] == 'A'||O[t+2] =='a') ) {
                    alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log07;
                }

                /* adjustment based on probability distribution */
                start_freq=0;
                freq_id = 0;

                double sub_sum = 0;
                int sub_count = 0;

                if (t>=30) {
                    for(int i=-30; i<=30; i++) {
                        if (t+i+2 < len_seq)
                        {
                            start_freq -= train_ptr->stop1[cg][i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                        }
                    }
                } else {
                    for(int i=(-1*t); i<=30; i++) {
                        if (t+i+2 < len_seq)
                        {
                            sub_sum += train_ptr->stop1[cg][i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                        }
                    }
                    sub_sum = sub_sum * 61 / (30 + t + 1);
                    start_freq -= sub_sum;
                }

                h_kd = train_ptr->E1_dist[cg][2] * exp(-1*pow(start_freq-train_ptr->E1_dist[cg][1],2)/(2*pow(train_ptr->E1_dist[cg][0],2)));
                r_kd = train_ptr->E1_dist[cg][5] * exp(-1*pow(start_freq-train_ptr->E1_dist[cg][4],2)/(2*pow(train_ptr->E1_dist[cg][3],2)));
                p_kd = h_kd / (h_kd + r_kd);

                if (p_kd<0.01) {
                    p_kd=0.01;
                } else if (p_kd>0.99) {
                    p_kd=0.99;
                }
                alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(p_kd);
            }
        }
        if (num_N>9) {
            for (int i=0; i<NUM_STATE; i++) if (i!=R_STATE) {
                    alpha[i][t] = max_dbl;
                    path[i][t] = R_STATE;
                }
        }
    }

    /***********************************************************/
    /* backtrack array to find the optimal path                */
    /***********************************************************/

    fprintf(fp_out, "%s\n", head);

    /* find the state for O[N] with the highest probability */
    for (int i = 0; i < NUM_STATE; i++) if (alpha[i][len_seq-1] < max_dbl) {
            max_dbl = alpha[i][len_seq-1];
            vpath[len_seq-1] = i;
        }

    /* backtrack the optimal path */
    for(int t=len_seq-2; t>=0; t--) vpath[t] = path[vpath[t+1]][t+1];

    int codon_start=0;
    start_t=-1;

    char codon[4], utr[65];
    char *dna = malloc(300000*sizeof(char));
    char *dna1 = malloc(300000*sizeof(char));
    char *dna_f = malloc(300000*sizeof(char));
    char *protein = malloc(100000*sizeof(char));
    int dna_id=0,dna_f_id=0,frame;
    for (int t=0; t<len_seq; t++) {

        if (codon_start==0 && start_t < 0 && ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) || (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1) || vpath[t] == S_STATE || vpath[t] == S_STATE_1 )) {
            dna_start_t_withstop=dna_start_t=start_t=t+1;
        }

        if (codon_start==0 && (vpath[t]==M1_STATE || vpath[t]==M4_STATE || vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)) {

            insert_id = 0;
            delete_id = 0;
            dna_id = 0;
            dna_f_id = 0;
            dna[dna_id]=O[t];
            dna_start_t_withstop=dna_start_t = t + 1; //Ye April 21, 2016
            if(t > 2 && (vpath[t] == M1_STATE_1 || vpath[t] == M4_STATE_1)) dna_start_t_withstop = t - 2;
            dna_f[dna_f_id]=O[t];
            start_orf=t+1;
            prev_match = vpath[t];

            codon_start = (vpath[t] < M6_STATE)?1:-1;

        } else if (codon_start!=0 && (vpath[t]==E_STATE || vpath[t]==E_STATE_1 || t==len_seq-1)) {

            if (vpath[t]==E_STATE || vpath[t]==E_STATE_1) {
                end_t=t+3;
            } else {
                /* FGS1.12 start: remove incomplete codon */
                temp_t = t;
                while(vpath[temp_t] != M1_STATE && vpath[temp_t] != M4_STATE  && vpath[temp_t] != M1_STATE_1  && vpath[temp_t] != M4_STATE_1) {
                    dna_f[dna_f_id] = '\0';
                    dna_f_id--;

                    dna[dna_id] = '\0';
                    dna_id--;

                    temp_t--;
                }
                end_t = temp_t;
                /* FGS1.12 end: remove incomplete codon */
            }

            if (dna_id > gene_len  ) {
                //these three lines moved here from outside of the loop above, YY July 23, 2018
                final_score = (alpha[vpath[end_t-4]][end_t-4]- alpha[vpath[start_t+2]][start_t+2] )/(end_t-start_t-5);
                frame = start_orf%3;
                if (frame==0) frame=3;

                if (codon_start==1) {
                    if(start_t == dna_start_t - 3) { //add complete start codon to dna, Ye April 21, 2016
                        dna_start_t -= 3;
                    }
                    if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
                        int start_old = start_t;
                        codon[0] = 0;
                        strncpy(codon, O + start_old-1, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save = 0;
                        int s_save = 0; //initialization, YY July 25 2018
                        while((!(!strcmp(codon, "TAA") || !strcmp(codon, "TAG") || !strcmp(codon, "TGA"))) && (start_old-1-s-35>=0)) {
                            if(!strcmp(codon, "ATG") || !strcmp(codon, "GTG") || !strcmp(codon, "TTG")) {
                                utr[0] = 0;
                                strncpy(utr, O+start_old-1-s-30,63);
                                utr[63] = 0;
                                //printf("check s=%d, codon %s\n", s, codon);
                                double freq_sum = 0;
                                for(int j = 0; j < strlen(utr) - 2; j ++) {
                                    int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]);
                                    freq_sum -= train_ptr->start[cg][j][idx];
                                    //printf("j=%d, key=%c%c%c %d, start %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->start[cg][j][idx]);
                                }
                                if(s == 0) {
                                    e_save = freq_sum;
                                    s_save = 0;
                                }
                                else if(freq_sum < e_save) {
                                    e_save = freq_sum;    //positive chain, upstream s_save = -1 * 3
                                    s_save = -1 * s;
                                }
                                //printf("s=%d freq_sum %lf\n", s, freq_sum);
                                //getchar();
                            }
                            s += 3;
                            codon[0] = 0;
                            strncpy(codon, O+start_old-1-s, 3);
                            codon[3] = 0;
                        }
                        //update start_t YY July 2018
                        if(s_save != 0) {
                            start_t = start_old+s_save;
                            dna_start_t = dna_start_t+s_save;
                        }
                    }

                    dna_end_t = end_t;
                    fprintf(fp_out, "%d\t%d\t+\t%d\t%lf\t", dna_start_t, dna_end_t, frame, final_score);
                    fprintf(fp_out, "I:");
                    for (int i=0; i<insert_id; i++) {
                        fprintf(fp_out, "%d,", insert[i]);
                    }
                    fprintf(fp_out, "\tD:");
                    for (int i=0; i<delete_id; i++) {
                        fprintf(fp_out, "%d,", delete[i]);
                    }
                    fprintf(fp_out, "\n");

                    //update dna before calling get_protein, YY July 2018
                    dna[0] = '\0';
                    strncpy(dna, O + dna_start_t - 1, dna_end_t - dna_start_t + 1);
                    dna[dna_end_t - dna_start_t + 1] = '\0';
                    //end of update dna
                    
                    get_protein(dna,protein,1, whole_genome);

                    fprintf(fp_aa, "%s_%d_%d_+\n", head, dna_start_t, dna_end_t);
                    fprintf(fp_dna, "%s_%d_%d_+\n", head, dna_start_t, dna_end_t);
                    fprintf(fp_aa, "%s\n", protein);
                    if (format==0) {
                        fprintf(fp_dna, "%s\n", dna);
                    } else if (format==1) {
                        fprintf(fp_dna, "%s\n", dna_f);
                    }
                } else if (codon_start==-1) {
                    if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
                        int end_old = end_t; //reverse
                        codon[0] = 0;
                        strncpy(codon, O + end_t-1-2, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save = 0;
                        int s_save = 0; //initialization, YY July 25, 2018
                        while((!(!strcmp(codon, "TTA") || !strcmp(codon, "CTA") || !strcmp(codon, "TCA"))) && (end_old-2+s+35 < len_seq) ) {
                            if(!strcmp(codon, "CAT") || !strcmp(codon, "CAC") || !strcmp(codon, "CAA")) {
                                utr[0] = 0;
                                strncpy(utr, O+end_old-1-2+s-30,63);
                                utr[63] = 0;
                                //printf("check s=%d, codon %s\n", s, codon);
                                double freq_sum = 0;
                                for(int j = 0; j < strlen(utr) - 2; j ++) {
                                    int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]);
                                    freq_sum -= train_ptr->stop1[cg][j][idx]; //stop1?? Ye, April 18, 2016
                                    //printf("j=%d, key=%c%c%c %d, stop1 %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->stop1[cg][j][idx]);
                                }
                                if(s == 0) {
                                    e_save = freq_sum;
                                    s_save = s;
                                }
                                else if(freq_sum < e_save) {
                                    e_save = freq_sum;    //negative chain, s_save = s, add, YY July 2018
                                    s_save = s;
                                }
                            }
                            s += 3;
                            codon[0] = 0;
                            strncpy(codon, O+end_old-1-2+s, 3);
                            codon[3] = 0;
                        }
                        //update end_t
                        end_t = end_old+s_save;
                    }

                    dna_end_t = end_t;
                    fprintf(fp_out, "%d\t%d\t-\t%d\t%lf\t", dna_start_t_withstop, dna_end_t, frame, final_score);
                    fprintf(fp_out, "I:");
                    for (int i=0; i<insert_id; i++) fprintf(fp_out, "%d,", insert[i]);
                    fprintf(fp_out, "\tD:");
                    for (int i=0; i<delete_id; i++) fprintf(fp_out, "%d,", delete[i]);
                    fprintf(fp_out, "\n");

                    //update dna before calling get_protein, YY July 2018
                    //use dna_end_t & dna_start_w_withstop to avoid incomplete codons & include start/stop codons
                    dna[0] = '\0';
                    strncpy(dna, O + dna_start_t_withstop - 1, dna_end_t - dna_start_t_withstop + 1);
                    dna[dna_end_t - dna_start_t_withstop + 1] = '\0';
                    //end of update dna
                    
                    
                    get_protein(dna,protein,-1, whole_genome); 

                    fprintf(fp_aa, "%s_%d_%d_-\n", head, dna_start_t_withstop, dna_end_t);
                    fprintf(fp_dna, "%s_%d_%d_-\n", head, dna_start_t_withstop, dna_end_t);

                    get_rc_dna(dna, dna1);
                    
                    fprintf(fp_aa, "%s\n", protein);
                    if (format==1) {
                        char *dna_f1 = malloc(300000*sizeof(char));
                        get_rc_dna_indel(dna_f, dna_f1);
                        fprintf(fp_dna, "%s\n", dna_f1);
                        free(dna_f1);
                    } else {
                        fprintf(fp_dna, "%s\n", dna1);
                    }
                }
            }
            codon_start=0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        } else if (codon_start!=0 &&
                   ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
                    (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1)) &&
                   vpath[t]-prev_match<6) {

            int out_nt = (vpath[t] < prev_match)? vpath[t]+6-prev_match:vpath[t]-prev_match;

            for (int kk=0; kk<out_nt; kk++) {  /* for deleted nt in reads */
                dna_id ++;
                dna[dna_id] = 'N';
                //printf("dna_id %d, dna-len %d\n", dna_id, strlen(dna));
                dna_f_id ++;
                dna_f[dna_f_id] = 'x';
                if (kk>0) {
                    delete[delete_id]=t+1;
                    delete_id++;
                }
            }
            dna[dna_id]=O[t];
            //printf("dna_id %d, add %d %c dna-len %d\n", dna_id, t, O[t], strlen(dna));
            dna_f[dna_f_id]=O[t];
            prev_match = vpath[t];

        } else if (codon_start!=0 &&
                   ((vpath[t]>=I1_STATE && vpath[t]<=I6_STATE) ||
                    (vpath[t]>=I1_STATE_1 && vpath[t]<=I6_STATE_1))) {
            dna_f_id ++;
            dna_f[dna_f_id] = tolower(O[t]);
            insert[insert_id]=t+1;
            insert_id++;

        }
        else if (codon_start!=0 && vpath[t]==R_STATE) {
            /* for long NNNNNNNNN, pretend R state */
            codon_start=0;
            start_t=-1;
            end_t = -1;
            dna_id=0;
            dna_f_id=0;

        }
    }

    free_dmatrix(alpha, NUM_STATE);
    free_imatrix(path, NUM_STATE);
    free_ivector(vpath);
    free(dna);
    free(dna1);
    free(dna_f);
    memset(protein,0, 100000);
    free(protein);
}

/* get_prob_from_cg
* This function determines the cg-ratio in a sequence.
* To determine the best chances of some transitions.
*
*/
int get_prob_from_cg(HMM *hmm_ptr, TRAIN *train_ptr, char *O) {
    int cg=0, len_seq= strlen(O);
    for (int i=0; i<len_seq; i++) {
        if ((O[i] == 'C'||O[i] =='c') || (O[i] == 'G'||O[i] == 'g') ) cg++;
    }
    cg = floor((cg*1.0/len_seq)*100)-26;
    return (cg < 0)?0:((cg > 43)?43:cg);
}


/**
* Reads files.
* 1. Reads trasition file and store in hmm datastructure
*/
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename,
                         char *sfilename,char *pfilename,char *s1filename,char *p1filename,char *dfilename, TRAIN *train_ptr) {

    double prob;
    char name[10];
    char start[10];
    char end[10];

    /* probabilities saved in log Ye April 18, 2016 */
    /* start <- ./train/start (start in forward)
       stop <- ./train/stop (stop in forward)
       start1 <- ./train/stop1 (start in reverse)
       stop1 <- ./train/start1 (stop in reverse)
    */

    /*
    * read all states, filename = train file */
    /****************************************************/
    /* transition                                       */
    /****************************************************/
    FILE *fp = fopen (filename, "r");

    /* Transition */
    fscanf(fp, "%*s");
    for (int i=0; i<14; i++) {
        fscanf(fp, "%s %lf", name, &prob);
        hmm_ptr->tr[tr2int(name)] = log(prob);
    }

    /* TransitionMI */
    fscanf(fp, "%*s");
    for (int i=0; i<16; i++) {
        fscanf(fp, "%s %s %lf\n", start, end, &prob);
        hmm_ptr->tr_M_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
    }

    /* TransitionII */
    fscanf(fp, "%*s");
    for (int i=0; i<16; i++) {
        fscanf(fp, "%s %s %lf", start, end, &prob);
        hmm_ptr->tr_I_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
    }

    /* PI */
    fscanf(fp, "%*s");
    for (int i=0; i<NUM_STATE; i++) {
        fscanf(fp, "%s %lf", name, &prob);
        hmm_ptr->pi[i] = log(prob);
    }
    fclose(fp);


    /****************************************************/
    /* M state transition                               */
    /****************************************************/
    FILE *fpm = fopen (mfilename, "r");
    for (int p=0; p<44; p++) {                       /* cg */
        fscanf(fp, "%*s");
        for (int i=0; i<6; i++) {                      /* period */
            for (int j=0; j<16; j++) {                   /* condition */
                for (int k=0; k<4; k++) {                  /* emission */
                    fscanf(fpm, "%lf", &prob);
                    train_ptr->trans[p][i][j][k] = log(prob);
                }
            }
        }
    }
    fclose(fpm);


    /****************************************************/
    /* M state_1 transition                             */
    /****************************************************/
    FILE * fpm1 = fopen (mfilename1, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int i=0; i<6; i++) {
            for (int j=0; j<16; j++) {
                for (int k=0; k<4; k++) {
                    fscanf(fpm1, "%lf", &prob);
                    train_ptr->rtrans[p][i][j][k] = log(prob);
                }
            }
        }
    }
    fclose(fpm1);


    /****************************************************/
    /* noncoding state  transition                      */
    /****************************************************/
    FILE *fpn = fopen (nfilename, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int j=0; j<4; j++) {
            for (int k=0; k<4; k++) {
                fscanf(fpn, "%lf", &prob);
                train_ptr->noncoding[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpn);


    /****************************************************/
    /* start                                            */
    /****************************************************/
    FILE * fps = fopen (sfilename, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int j=0; j<61; j++) {
            for (int k=0; k<64; k++) {
                fscanf(fps, "%lf", &prob);
                train_ptr->start[p][j][k] = log(prob);
            }
        }
    }
    fclose(fps);


    /****************************************************/
    /* stop                                             */
    /****************************************************/
    FILE *fpp = fopen (pfilename, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int j=0; j<61; j++) {
            for (int k=0; k<64; k++) {
                fscanf(fpp, "%lf", &prob);
                train_ptr->stop[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpp);


    /****************************************************/
    /* start1                                           */
    /****************************************************/
    FILE *fps1 = fopen (s1filename, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int j=0; j<61; j++) { //58->61 Ye, April 18, 2016
            for (int k=0; k<64; k++) {
                fscanf(fps1, "%lf", &prob);
                train_ptr->start1[p][j][k] = log(prob);
            }
        }
    }
    fclose(fps1);


    /****************************************************/
    /* stop1                                            */
    /****************************************************/
    FILE *fpp1 = fopen (p1filename, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int j=0; j<61; j++) {
            for (int k=0; k<64; k++) {
                fscanf(fpp1, "%lf", &prob);
                train_ptr->stop1[p][j][k] = log(prob);
            }
        }
    }
    fclose(fpp1);


    /****************************************************/
    /* pwm distribution                                 */
    /* S_dist, E_dist, S1_dist, E1_dist NOT in log      */
    /****************************************************/
    FILE *fpd = fopen (dfilename, "r");
    for (int p=0; p<44; p++) {
        fscanf(fp, "%*s");
        for (int k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->S_dist[p][k] = prob;
        }
        for (int k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->E_dist[p][k] = prob;
        }
        for (int k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->S1_dist[p][k] = prob;
        }
        for (int k=0; k<6; k++) {
            fscanf(fpd, "%lf", &prob);
            train_ptr->E1_dist[p][k] = prob;
        }
    }
    fclose(fpd);

}