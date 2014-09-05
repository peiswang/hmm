/*
 **      Author: peisong wang peisong.wang@nlpr.ia.ac.cn
 **      Date:   5 September 2014 
 **      File:   forward.c
 **      Purpose: Forward algorithm for computing the probabilty
 **              of observing a sequence given a HMM model parameter.
 **      Organization: 
 **
 **      AN Extension of UMDHMM
 */
#include <stdio.h>
#include <float.h>
#include "nrutil.h"
#include "hmm.h"


void Forward(HMM *phmm, int T, double **outprob, double **alpha, double *pprob)
{
    int     i, j;   /* state indices */
    int     t;      /* time index */

    double sum;     /* partial sum */
 
    /* 1. Initialization */
 
    for (i = 1; i <= phmm->N; i++)
    {
        alpha[1][i] = (phmm->pi[i] + outprob[i][1]);
    }
 
    /* 2. Induction */
 
    for (t = 1; t < T; t++) {
        for (j = 1; j <= phmm->N; j++) {
            sum = NAN;
            for (i = 1; i <= phmm->N; i++)
                sum = logSum( sum, (alpha[t][i] + phmm->A[i][j]) );

            alpha[t+1][j] = (sum + outprob[j][t+1]);
        }
    }
 
    /* 3. Termination */
    *pprob = NAN;
    for (i = 1; i <= phmm->N; i++)
        *pprob = logSum(*pprob, alpha[T][i]);
 
}

void ForwardWithScale(HMM *phmm, int T, double **outprob, double **alpha, 
	double *scale, double *pprob)
/*  pprob is the LOG probability */
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
    double c_s; /* 1/scale */

	/* 1. Initialization */

	scale[1] = 0.0;	
	for (i = 1; i <= phmm->N; i++) {
		alpha[1][i] = phmm->pi[i] * outprob[i][1];
		scale[1] += alpha[1][i];
	}
    if(isnan(scale[1]) || scale[1]<DBL_MIN)
    {
        printf("Scale[1] is 0 in ForwardWithScale()!!!\n");
        return;
    }
    c_s = 1.0 / scale[1];

	for (i = 1; i <= phmm->N; i++) 
		alpha[1][i] *= c_s;
	
	/* 2. Induction */

	for (t = 1; t <= T - 1; t++) {
		scale[t+1] = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			sum = 0.0;
			for (i = 1; i <= phmm->N; i++) 
				sum += alpha[t][i]* (phmm->A[i][j]); 

			alpha[t+1][j] = sum * outprob[j][t+1];
			scale[t+1] += alpha[t+1][j];
		}
        if(scale[t+1]<DBL_MIN)
        {
            printf("Scale[%d] is 0 in ForwardWithScale()!!!\n", t+1);
            return;
        }
        c_s = 1.0 / scale[t+1];

		for (j = 1; j <= phmm->N; j++) 
			alpha[t+1][j] *= c_s; 
	}

	/* 3. Termination */
	*pprob = 0.0;

	for (t = 1; t <= T; t++)
		*pprob += log(scale[t]);
	
}
