/*
**      Author: peisong wang peisong.wang@nlpr.ia.ac.cn
**      Date:   5 September 2014 
**      File:   backward.c
**      Purpose: Backward algorithm for computing the probabilty
**              of observing a sequence given a HMM model parameter.
**      Organization: 
**
**      AN Extension of UMDHMM
*/

#include <stdio.h>
#include "hmm.h"
#include "nrutil.h"

void Backward(HMM *phmm, int T, double **outprob, double **beta, double *pprob)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
        double sum;
 
 
        /* 1. Initialization */
 
        for (i = 1; i <= phmm->N; i++)
                beta[T][i] = 0.0;
 
        /* 2. Induction */
 
        for (t = T - 1; t >= 1; t--) {
                for (i = 1; i <= phmm->N; i++) {
                        sum = NAN;
                        for (j = 1; j <= phmm->N; j++)
                                sum = logSum(sum, phmm->A[i][j]+outprob[j][t+1]+beta[t+1][j]);
                        beta[t][i] = sum;
                }
        }
 
        /* 3. Termination */
        *pprob = NAN;
        for (i = 1; i <= phmm->N; i++)
                *pprob = logSum(*pprob, phmm->pi[i]+outprob[i][1]+beta[1][i]);
 
}

void BackwardWithScale(HMM *phmm, int T, double **outprob, double **beta, 
	double *scale, double *pprob)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
        /* 1. Initialization */
 
        for (i = 1; i <= phmm->N; i++)
                beta[T][i] = 1.0;//scale[T]; 
 
        /* 2. Induction */
 
        for (t = T - 1; t >= 1; t--) {
            for (i = 1; i <= phmm->N; i++) 
            {
			    sum = 0.0;
                for (j = 1; j <= phmm->N; j++)
                    sum += phmm->A[i][j] * outprob[j][t+1]*beta[t+1][j];
                beta[t][i] = sum/scale[t+1];
 
            }
        }
 
}
