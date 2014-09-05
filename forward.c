/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   forward.c
**      Purpose: Foward algorithm for computing the probabilty 
**		of observing a sequence given a HMM model parameter.
**      Organization: University of Maryland
**
**      $Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $
*/
#include <stdio.h>
#include <float.h>
#include "nrutil.h"
#include "hmm.h"

static char rcsid[] = "$Id: forward.c,v 1.2 1998/02/19 12:42:31 kanungo Exp kanungo $";

//inline double f(double x, double y)
//{
//    return x+y;
//}

void Forward(HMM *phmm, int T, double **outprob, double **alpha, double *pprob)
{
    int     i, j;   /* state indices */
    int     t;      /* time index */

    double sum;     /* partial sum */
 
    /* 1. Initialization */
 
    for (i = 1; i <= phmm->N; i++)
    {
        //printf("pi[%d] = %.5f, outprob[%d][1] = %.5f\n ", i, phmm->pi[i],i,outprob[i][1]);
        alpha[1][i] = (phmm->pi[i] + outprob[i][1]);
        //printf("%.5f  \n",logProd(phmm->pi[i], outprob[i][1]));
        //printf("%.5f  ",alpha[1][i]);
        //printf("%.5f  \n",f(phmm->pi[i],outprob[i][1]));
    }
    //printf("\n");
 
    /* 2. Induction */
 
    for (t = 1; t < T; t++) {
        for (j = 1; j <= phmm->N; j++) {
            sum = NAN;
            for (i = 1; i <= phmm->N; i++)
                sum = logSum( sum, (alpha[t][i] + phmm->A[i][j]) );

            alpha[t+1][j] = (sum + outprob[j][t+1]);
            //printf("%.5f  ",alpha[t+1][j]);
        }
    //printf("\n");
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
    if(isnan(scale[1]))
        printf("sssssssssssssssss\n");
    if(scale[1]<DBL_MIN)
    {
        printf("Scale[1] is 0 in ForwardWithScale()!!!\n");
        //fprintf(stderr, "Scale[1] is 0 in ForwardWithScale()!!!\n");
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
            //fprintf(stderr, "Scale[%d] is 0 in ForwardWithScale()!!!\n", t+1);
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
