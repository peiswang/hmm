/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   baumwelch.c
**      Purpose: Baum-Welch algorithm for estimating the parameters
**              of a HMM model, given an observation sequence. 
**      Organization: University of Maryland
**
**	Update: 
**	Author: Tapas Kanungo
**	Date:	19 April 1999
**	Purpose: Changed the convergence criterion from ratio
**		to absolute value. 
**
**      $Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $
*/

#include <stdio.h> 
#include <float.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"
#include "sample.h"

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA -100 

struct local_store_c
{
    double **alpha;
    double **beta;
    double *scale;
    double **outprob;
};

struct local_store_a
{
    double **numeratorA;
    double *denominatorA;
    double **numeratorMiu;
    double **numeratorCov;
    double *denominatorM;
};

void f(HMM *phmm, struct samples *p_samples, struct local_store_c *c, struct local_store_a *acc, double *logp)
{	
    int i, j, k, t;
    int sample_index;

    int T;
    int N = phmm->N;
    int D = phmm->D;
    int sample_count = p_samples->sample_count;
    int *fps = p_samples->feature_count_per_sample;

    double **alpha = c->alpha;
    double **beta = c->beta;
    double *scale = c->scale;
    double **outprob = c->outprob;

    double **sample;
    double logprobf, logprobb;
    *logp = NAN;
    *logp = 0;
    int valid_count = 0;

    for(sample_index=0;sample_index<sample_count;sample_index++)
    {
        T = fps[sample_index];
        sample = p_samples->data[sample_index];
        CalcOutProb(phmm, sample, T, outprob);
	    Forward(phmm, T, outprob, alpha, &logprobf);
	    Backward(phmm, T, outprob, beta, &logprobb);

        //if(logprobf<600)
        //        continue;
        valid_count++;
        printf("logprobf is %lf\n", logprobf);

        if(fabs(logprobf-logprobb) > 0.00000001)
                printf("logprobf != logrpbb \n");
        //printf("logprobf = %f\n",(logprobf));

        //*logp = logSum(*logp, logprobf);
        *logp += logprobf;

        /* update pi[]*/
        //for(i=1;i<N;i++)
        //    phmm->pi[i] = 0.0;
        //for(sample_id=1;sample_id<=sample_count;sample_id++)
        //{
        //    for(i=1;i<=N;i++)
        //        phmm->pi[i] += .001 + .999 * gamma[sample_id][1][i];
        //}
        //for(i=1;i<N;i++)
        //    phmm->pi[i] /= sample_count;

        /* accumulate A*/
        for(t=1;t<=T-1;t++)
        {
            //double scale_inv = 1.0 / scale[t+1];
            for(i=1;i<=N;i++)
            {
                for(j=1;j<=N;j++)
                {
                    if(phmm->A_topo[i][j]==0)
                        continue;
                    double xi = ( alpha[t][i] + phmm->A[i][j] + 
                                         outprob[j][t+1] + beta[t+1][j]-logprobf
                                );
                    acc->numeratorA[i][j] = logSum(acc->numeratorA[i][j], xi);
                    acc->denominatorA[i] = logSum(acc->denominatorA[i], xi);
                }
            }
        }

        for(j=1;j<=N;j++)
        {
            for(t=1;t<=T;t++)
            {
                double gamma_t_j = NAN;
                if(t==1)
                    gamma_t_j = phmm->pi[j];
                else
                {
                    for(i=1;i<=N;i++)
                    {
                        gamma_t_j = logSum(gamma_t_j, (alpha[t-1][i] + phmm->A[i][j]));
                    }
                }
                //printf("gamma is %f\n", scale[t]);
                gamma_t_j = (gamma_t_j + outprob[j][t] + beta[t][j] - logprobf);

                gamma_t_j = logExp(gamma_t_j);
                //printf("gamma_t_j is %.20f\n",gamma_t_j);

                for(k=1;k<=D;k++)
                {
                    acc->numeratorMiu[j][k] += gamma_t_j * sample[t-1][k-1]; 

                    acc->numeratorCov[j][k] += gamma_t_j * (sample[t-1][k-1]-phmm->B.miu[j][k]) *
                                                              (sample[t-1][k-1]-phmm->B.miu[j][k]);
                }
                acc->denominatorM[j] += gamma_t_j;
            }
        }
    } /* for sample*/

    for(i=1;i<=N;i++)
    {
        /* update A*/
        if(i<N)
        {
            ////////////////  denoninator < -10000000
            if(isnan(acc->denominatorA[i]))
            {
                printf("denominatorA[%d] is 0!!!\n", i);
            }

            for(j=1;j<=N;j++)
            {
                phmm->A[i][j] = acc->numeratorA[i][j] - acc->denominatorA[i];
                if(isnan(phmm->A[i][j]) && phmm->A_topo[i][j]==1)
                {
                    printf("A[%d][%d] is nan\n", i, j);
                }
                if(phmm->A[i][j]>0)
                {
                        
                    printf("A[%d][%d] = %f - %f =  %f\n", i, j, acc->numeratorA[i][j],
                                    acc->denominatorA[i],(phmm->A[i][j]));
                    phmm->A[i][j] = 0.0;
                }
            }
        }

        /* update B*/
        
        double c_d;

        if(acc->denominatorM[i]<DBL_MIN)
        {
            printf("denominatorM[%d] is 0!!!\n", i);
        }
        else
            c_d = 1.0 / acc->denominatorM[i];

        for(k=1;k<=D;k++)
        {
            phmm->B.miu[i][k] = acc->numeratorMiu[i][k] * c_d;
            phmm->B.cov[i][k] = acc->numeratorCov[i][k] * c_d + MIN_COV;
            phmm->B.cov_inv[i][k] = 1.0 / phmm->B.cov[i][k];
        }

    }

    //*logp /= sample_count;
    *logp /= valid_count;
}

void BaumWelch(HMM *phmm, struct samples *p_samples, int *piter, 
	double *plogprobinit, double *plogprobfinal, int maxiter)
{
    int D = phmm->D;
    int N = phmm->N;
    int max_t = p_samples->feature_count_max;
    int sample_count = p_samples->sample_count;

    struct local_store_c cs;
    struct local_store_a as;

    double **outprob;
	double delta, logprob, logprobprev;

	//double deltaprev = 10e-70;

    cs.scale = dvector(1, max_t);
    cs.alpha = dmatrix(1, max_t, 1, N);
    cs.beta = dmatrix(1, max_t, 1, N);
    cs.outprob = dmatrix(1, phmm->N, 1, max_t);

    as.numeratorA = dmatrix(1, N, 1, N);
    clear_dmatrix_nan(as.numeratorA, 1, N, 1, N);
    as.denominatorA = dvector(1, N);
    clear_dvector_nan(as.denominatorA, 1, N);
    as.numeratorMiu = dmatrix(1, N, 1, D);
    clear_dmatrix(as.numeratorMiu, 1, N, 1, D);
    as.numeratorCov = dmatrix(1, N, 1, D);
    clear_dmatrix(as.numeratorCov, 1, N, 1, D);
    as.denominatorM = dvector(1, N);
    clear_dvector(as.denominatorM, 1, N);


    /* calc output probability after miu and cov are changed*/
    logprobprev = -1000;
    
    *piter = 0;
    //while( *piter < 1)
    while( *piter < maxiter )
    {
        *piter = *piter + 1;
	    /* compute difference between log probability of 
	          two iterations */
        f(phmm, p_samples, &cs, &as, &logprob);

	    delta = logprob - logprobprev; 
	    logprobprev = logprob;
        printf("iter %d, delta is : %.20f\n", *piter, delta);
        //printf("iter %d, prob is : %.20f\n", *piter, logprob);
        if(delta<DELTA)
            break;
    }
	//while (delta > DELTA); 
 
	*plogprobfinal = logprob; /* log P(O|estimated model) */
	free_dvector(cs.scale, 1, max_t);
    free_dmatrix(cs.alpha, 1, max_t, 1, N);
    free_dmatrix(cs.beta, 1, max_t, 1, N);
    free_dmatrix(cs.outprob, 1, phmm->N, 1, max_t);
    free_dmatrix(as.numeratorA, 1, N, 1, N);
	free_dvector(as.denominatorA, 1, N);
    free_dmatrix(as.numeratorMiu, 1, N, 1, D);
    free_dmatrix(as.numeratorCov, 1, N, 1, D);
	free_dvector(as.denominatorM, 1, N);


}

