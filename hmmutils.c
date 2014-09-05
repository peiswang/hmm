/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   hmmutils.c
**      Purpose: utilities for reading, writing HMM stuff. 
**      Organization: University of Maryland
**
**      $Id: hmmutils.c,v 1.4 1998/02/23 07:51:26 kanungo Exp kanungo $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "nrutil.h"
#include "hmm.h"
#include "sample.h"

static char rcsid[] = "$Id: hmmutils.c,v 1.4 1998/02/23 07:51:26 kanungo Exp kanungo $";

void ReadHMM(FILE *fp, HMM *phmm)
{
	int i, j, k, l;

	fscanf(fp, "N= %d\n", &(phmm->N)); 

	fscanf(fp, "M= %d\n", &(phmm->M)); 
    
	fscanf(fp, "D= %d\n", &(phmm->D)); 

	fscanf(fp, "A:\n");
	phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++) { 
		for (j = 1; j <= phmm->N; j++) {
			fscanf(fp, "%lf", &(phmm->A[i][j])); 
		}
		fscanf(fp,"\n");
	}

	fscanf(fp, "A_topo:\n");
	phmm->A_topo = (int **) imatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++) { 
		for (j = 1; j <= phmm->N; j++) {
			fscanf(fp, "%d", &(phmm->A_topo[i][j])); 
		}
		fscanf(fp,"\n");
	}

	fscanf(fp, "B_miu:\n");
	phmm->B.miu = (double **) dmatrix(1, phmm->N, 1, phmm->D);
	for (j = 1; j <= phmm->N; j++) { 
		for (k = 1; k <= phmm->D; k++) {
			fscanf(fp, "%lf", &(phmm->B.miu[j][k])); 
		}
		fscanf(fp,"\n");
	}

	fscanf(fp, "B_cov:\n");
	phmm->B.cov= (double **) dmatrix(1, phmm->N, 1, phmm->D);
	for (j = 1; j <= phmm->N; j++) { 
		for (k = 1; k <= phmm->D; k++) {
			fscanf(fp, "%lf", &(phmm->B.cov[j][k])); 
		}
		fscanf(fp,"\n");
	}

	fscanf(fp, "pi:\n");
	phmm->pi = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++) 
		fscanf(fp, "%lf", &(phmm->pi[i])); 

    /* allocate cov_inv*/
	phmm->B.cov_inv = (double **) dmatrix(1, phmm->N, 1, phmm->D);
	for (j = 1; j <= phmm->N; j++) { 
		for (k = 1; k <= phmm->D; k++) {
			phmm->B.cov_inv[j][k] = 1.0 / phmm->B.cov[j][k]; 
		}
	}

}

void FreeHMM(HMM *phmm)
{
	free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_imatrix(phmm->A_topo, 1, phmm->N, 1, phmm->N);
	free_dmatrix(phmm->B.miu, 1, phmm->N, 1, phmm->D);
	free_dmatrix(phmm->B.cov, 1, phmm->N, 1, phmm->D);
	free_dmatrix(phmm->B.cov_inv, 1, phmm->N, 1, phmm->D);
	free_dvector(phmm->pi, 1, phmm->N);

    /* wrong*/
	//free_dmatrix(phmm->B.cov_inv, 1, phmm->N, 1, phmm->D);
	//free_dmatrix3d(phmm->B.outprob, 1, phmm-> 1, phmm->N, 1, phmm->D);
}

/*
** InitHMM() This function initializes matrices A, B and vector pi with
**	random values. Not doing so can result in the BaumWelch behaving
**	quite weirdly.
*/ 

void InitHMM(HMM *phmm, int N, int M, int D, struct samples *p_samples,int **topo)
{
    assert(N >= 1);
	int i, j, k;
	double sum;

    phmm->D = D;
    phmm->M = M;
    phmm->N = N;

    int *feature_count_per_sample = p_samples->feature_count_per_sample;
    int sample_count = p_samples->sample_count;
    double feature_count = 0.0;
    for(i=0;i<sample_count;i++)
    {
        feature_count += feature_count_per_sample[i];
    }

    double aveSampleLen = feature_count / sample_count;
    double stateTranProb = (N-1) / (aveSampleLen-1);
 
    /* initialize pi*/
	phmm->pi = (double *) dvector(1, phmm->N);
    for(i=1;i<=phmm->N;i++)
        phmm->pi[i] = NAN;
    //clear_dvector(phmm->pi, 1, phmm->N);
    phmm->pi[1] = 0.0;
    
    /* initialize A*/
    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    //clear_dmatrix(phmm->A, 1, N, 1, N);
    for(i=1;i<=phmm->N;i++)
    {
        for(j=1;j<=phmm->N;j++)
            phmm->A[i][j] = NAN;
    }
    for (i = 1; i < phmm->N; i++) // < N
    {
        phmm->A[i][i] = log(1 - stateTranProb);
        phmm->A[i][i+1] = log(stateTranProb);
	}
    phmm->A[phmm->N][phmm->N] = 0.0;

    /* initialize topo*/
    phmm->A_topo = (int **) imatrix(1, phmm->N, 1, phmm->N);
    clear_imatrix(phmm->A_topo, 1, N, 1, N);
    if(topo)
    {
        for (i = 1; i <= N; i++) 
        {
            for(j=1;j<=N;j++)
                phmm->A_topo[i][j], topo[i][j];
	    }
    }else
    {
        for (i = 1; i <= N; i++) 
        {
            for(j=1;j<=N;j++)
                if(j==i || j-i==1)
                    phmm->A_topo[i][j] = 1;
                else
                    phmm->A_topo[i][j] = 0;
	    }
    }
 
    /* initialize B*/
    phmm->B.miu = (double **) dmatrix(1, phmm->N, 1, phmm->D);
	phmm->B.cov = (double **) dmatrix(1, phmm->N, 1, phmm->D);
	phmm->B.cov_inv = (double **) dmatrix(1, phmm->N, 1, phmm->D);
    InitMC(phmm, p_samples);

    /* initialize some useful structure*/
    //outprob = (double ***) dmatrix3d(1, p_samples->sample_count, 1, phmm->N, 1, p_samples->feature_count_max);

}

void CopyHMM(HMM *phmm1, HMM *phmm2)
{
        int i, j, k;

        phmm2->N = phmm1->N;
 
        phmm2->M = phmm1->M;

        phmm2->D = phmm1->D;

        phmm2->A = (double **) dmatrix(1, phmm2->N, 1, phmm2->N);
 
        for (i = 1; i <= phmm2->N; i++)
                for (j = 1; j <= phmm2->N; j++)
                        phmm2->A[i][j] = phmm1->A[i][j];
 
        phmm2->B.miu = (double **) dmatrix(1, phmm2->N, 1, phmm2->D);
        for (j = 1; j <= phmm2->N; j++)
                for (k = 1; k <= phmm2->D; k++)
                        phmm2->B.miu[j][k] = phmm1->B.miu[j][k];
 
        phmm2->B.cov= (double **) dmatrix(1, phmm2->N, 1, phmm2->D);
        for (i = 1; i <= phmm2->N; i++)
                for (j = 1; j <= phmm2->D; j++)
                    phmm2->B.cov[i][j] = phmm1->B.cov[i][j];

        phmm2->B.cov_inv = (double **) dmatrix(1, phmm2->N, 1, phmm2->D);
        for (i = 1; i <= phmm2->N; i++)
                for (j = 1; j <= phmm2->D; j++)
                    phmm2->B.cov_inv[i][j] = phmm1->B.cov_inv[i][j];

        phmm2->pi = (double *) dvector(1, phmm2->N);
        for (i = 1; i <= phmm2->N; i++)
            phmm2->pi[i] = phmm1->pi[i]; 
 
}

void PrintHMM(const char *filename, HMM *phmm)
{
        int i, j, k;
    FILE *fp = fopen(filename,"w+");
    if(fp<0) 
        return;

	fprintf(fp, "N= %d\n", phmm->N); 
	fprintf(fp, "M= %d\n", phmm->M); 
	fprintf(fp, "D= %d\n", phmm->D); 
 
	fprintf(fp, "A:\n");
        for (i = 1; i <= phmm->N; i++) {
                for (j = 1; j <= phmm->N; j++) {
                        fprintf(fp, "%f ", logExp(phmm->A[i][j]) );
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "A_topo:\n");
        for (i = 1; i <= phmm->N; i++) {
                for (j = 1; j <= phmm->N; j++) {
                        fprintf(fp, "%d ", phmm->A_topo[i][j] );
		}
		fprintf(fp, "\n");
	}
 
	fprintf(fp, "B_miu:\n");
        for (j = 1; j <= phmm->N; j++) {
                for (k = 1; k <= phmm->D; k++){
                        fprintf(fp, "%f ", phmm->B.miu[j][k]);
		}
		fprintf(fp, "\n");
	}

    fprintf(fp, "B_cov:\n");
    for(i=1;i<=phmm->N;i++)
    {
        for(j=1;j<=phmm->D;j++)
        {
            fprintf(fp, "%f ", phmm->B.cov[i][j]);
        }
		fprintf(fp, "\n");
    }
 
	fprintf(fp, "pi:\n");
    for (i = 1; i <= phmm->N; i++) {
		fprintf(fp, "%f ", logExp(phmm->pi[i]));
	}
	fprintf(fp, "\n\n");
    fclose(fp);
}

void InitMC(HMM *phmm, struct samples * p_samples)
{
    int i, j, k;
    int N = phmm->N;
    int D = phmm->D;
    double **miu = phmm->B.miu;
    double **cov = phmm->B.cov;
    double **cov_inv = phmm->B.cov_inv;
    double *feature_tmp;
    clear_dmatrix(miu, 1, N, 1, D);
    clear_dmatrix(cov, 1, N, 1, D);

    int sample_count = p_samples->sample_count;
    int *fps = p_samples->feature_count_per_sample;
    int *index_count = (int*) ivector(1, N); //calloc( N, sizeof(int) );

    for(i=1;i<=N;i++)
        index_count[i] = 0;

    for(i=0;i<sample_count;i++)
    {
        for(j=0;j<fps[i];j++)
        {
            int index = round( j * N / fps[i]) + 1;
            index_count[index]++;
            feature_tmp = p_samples->data[i][j];
            for(k=1;k<=D;k++)
            {
                miu[index][k] += feature_tmp[k-1];
                cov[index][k] += feature_tmp[k-1] * feature_tmp[k-1];
            }
        
        }
    }

    for(i=1;i<=N;i++)
    {
        double c_c = 1.0 / index_count[i];
        for(j=1;j<=D;j++)
        {
            miu[i][j] *= c_c;
            cov[i][j] *= c_c;
            cov[i][j] -= miu[i][j] * miu[i][j];
            //if(cov[i][j]<MIN_COV)
            //    cov[i][j]=MIN_COV;
            cov[i][j] += MIN_COV;
            cov_inv[i][j] = 1.0 / cov[i][j];
        }
    }
    free_ivector(index_count, 1, N);
}

// calc outprobability after updating miu and cov
void CalcOutProb(HMM *phmm, double **sample, int T, double **outprob)
{
    int f, i, j;
    int N = phmm->N;
    int D = phmm->D;

    double *feature_tmp;

    double **miu = phmm->B.miu;
    double **cov = phmm->B.cov;
    double **cov_inv = phmm->B.cov_inv;
    //double **cov_inv = phmm->B.cov_inv;

    double prob1 = - 0.5 * N * log(2*M_PI); // -N*log(2pi)/2
    double *prob2 = (double*)dvector(1,N);  // -1/2*log(|sigma|)

    for(i=1;i<=N;i++)
    {
        double tmp = 0.0;
        for(j=1;j<=D;j++)
            tmp += log(cov[i][j]);
        prob2[i] = - tmp * 0.5;
        //printf("prob2 %f\n",prob2[i]);
    }
   

    for(f=0;f<T;f++)
    {
            double tmp, x;
            feature_tmp = sample[f];
            for(i=1;i<=N;i++)
            {
                tmp = 0.0;
                for(j=1;j<=D;j++)
                {
                    x = feature_tmp[j-1]-miu[i][j];
                    tmp += x*x*cov_inv[i][j];
                }
                //outprob[i][f+1] = exp(prob1 + prob2[i] -0.5 * tmp);
                outprob[i][f+1] = prob1 + prob2[i] - 0.5 * tmp;
                //printf("%.5f  ", prob1+prob2[i]-0.5*tmp);
                //printf("%.3f  ", outprob[i][f+1]);
            } 
            //printf("\n");
    }
    //printf("\n");
    free_dvector(prob2,1,N);

}

//void CalcOutProb(HMM *phmm, struct samples * p_samples, double ***outprob, double **cov_inv)
//{
//    int s, f, i, j;
//    int N = phmm->N;
//    int D = phmm->D;
//    int sample_count = p_samples->sample_count;
//    int *fps = p_samples->feature_count_per_sample;
//    double *feature_tmp;
//
//    double **miu = phmm->B.miu;
//    double **cov = phmm->B.cov;
//    //double **cov_inv = phmm->B.cov_inv;
//    for(i=1;i<=N;i++)
//        for(j=1;j<=D;j++)
//        {
//            if(cov[i][j]<DBL_MIN){
//                printf("cov[%d][%d] is 0!!!\n",i,j);
//            }
//            else{
//                cov_inv[i][j] = 1.0 / cov[i][j];
//            }
//        }
//
//
//    double prob1 = - 0.5 * N * log(2*M_PI); // -N*log(2pi)/2
//    double *prob2 = (double*)dvector(1,N);  // -1/2*log(|sigma|)
//    for(i=1;i<=N;i++)
//    {
//        double tmp = 0.0;
//        for(j=1;j<=D;j++)
//            tmp += log(cov[i][j]);
//        prob2[i] = - tmp * 0.5;
//    }
//
//    for(s=0;s<sample_count;s++)
//    {
//        for(f=0;f<fps[s];f++)
//        {
//            double tmp, x;
//            feature_tmp = p_samples->data[s][f];
//            for(i=1;i<=N;i++)
//            {
//                tmp = 0.0;
//                for(j=1;j<=D;j++)
//                {
//                    x = feature_tmp[j-1]-miu[i][j];
//                    tmp += x*x*cov_inv[i][j];
//                }
//                outprob[s+1][i][f+1] = exp(prob1 + prob2[i] -0.5 * tmp);
//                //printf("%.5f  ", prob1+prob2[i]-0.5*tmp);
//                //printf("%.5f  ", outprob[s+1][i][f+1]);
//            } 
//            //printf("\n");
//        }
//        //printf("\n");
//    }
//    printf("\n");
//    free_dvector(prob2,1,N);
//
//}
//
