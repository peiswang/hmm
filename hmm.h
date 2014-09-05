/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   hmm.h
**      Purpose: datastructures used for HMM. 
**      Organization: University of Maryland
**
**	Update:
**	Author: Tapas Kanungo
**	Purpose: include <math.h>. Not including this was
**		creating a problem with forward.c
**      $Id: hmm.h,v 1.9 1999/05/02 18:38:11 kanungo Exp kanungo $
*/

#ifndef __HMM_H__
#define __HMM_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sample.h"


#define MIN_COV     1E-4
//#define MIN_COV 0.00000001


typedef struct{
    double  **miu;    /* M means of mixture gaussian*/
    double  **cov;   /* M covariance of mixture gaussian*/

    double  **cov_inv;
    //double  ***outprob;
} output;

typedef struct {
	int N;		/* number of states;  Q={1,2,...,N} */
	int M; 		/* number of observation symbols; V={1,2,...,M}*/
    int D;    /* dimension of observation vector */
	double	**A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
    int **A_topo;
	output B;   /* B[1..N][1..M]. b[j][k] is the probability of
			   of observing symbol k in state j */
	double	*pi;	/* pi[1..N] pi[i] is the initial state distribution. */
} HMM;

void ReadHMM(FILE *fp, HMM *phmm);
void PrintHMM(const char *fn, HMM *phmm);
void AllocateHMM(HMM *phmm, int N, int M, int D);
void InitHMM(HMM *phmm, int N, int M, int D, struct samples *p_samples,int **topo);
void CopyHMM(HMM *phmm1, HMM *phmm2);
void FreeHMM(HMM *phmm);

//void ReadSequence(FILE *fp, int *pT, int **pO);
//void PrintSequence(FILE *fp, int T, int *O);
//void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q);
//int GenInitalState(HMM *phmm);
//int GenNextState(HMM *phmm, int q_t);
//int GenSymbol(HMM *phmm, int q_t);

void InitMC(HMM *phmm, struct samples * p_samples);
void CalcOutProb(HMM *phmm, double **sample, int T, double **outprob);

 
void Forward(HMM *phmm, int T, double **outprob, double **alpha, double *pprob);
void ForwardWithScale(HMM *phmm, int T, double **outprob, double **alpha,
        double *scale, double *pprob);
void Backward(HMM *phmm, int T, double **outprob, double **beta, double *pprob);
void BackwardWithScale(HMM *phmm, int T, double **outprob, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMM *phmm, struct samples *, 
           int *niter, double *plogprobinit, double *plogprobfinal, int maxiter);

//void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi,
//        int *q, double *pprob);
//void ViterbiLog(HMM *phmm, int T, int *O, double **delta, int **psi,
//        int *q, double *pprob);

/* random number generator related functions*/

//int hmmgetseed(void);
//void hmmsetseed(int seed);
//double hmmgetrand(void);

#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
 
#endif
