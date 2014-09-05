#include <stdio.h>
#include <math.h>
#include "hmm.h"
#include "nrutil.h"
#include "sample.h"
#include "rand.h"

HMM* createHMM(int N, int D);

int main()
{
        int i,j,k;
    HMM hmm;
    HMM *phmm = &hmm;



    int N = 5;
    int D = 1764;
    int T = 100;
    HMM *phmm_gt = createHMM(N,D);
    PrintHMM("test_gt.hmm",phmm_gt);
    FILE *ptest = fopen("a.data","w+");
    double **data = (double **) dmatrix(1,T,1,D);
    fprintf(ptest,"sample_count= 140\n"); 
    for(i=1;i<=140;i++)
    {
        //fprintf(ptest,"# %d\n",dd); 
        int state_count = 0;
        int last_state_count = 0;
        int current_state = 1;
        while((last_state_count==0)||state_count>N*last_state_count&&state_count<T)
        {
            state_count++;
            if(current_state==N)
                last_state_count++;

            for(k=1;k<=D;k++)
            {
                data[state_count][k] = getNormalRand(phmm_gt->B.miu[current_state][k],
                                sqrt(phmm_gt->B.cov[current_state][k]));
            }
            if(getRand()>logExp(phmm_gt->A[current_state][current_state]))
            {
                current_state++;
            }
        }
        fprintf(ptest,"# %d\n", state_count);
        for(j=1;j<=state_count;j++)
        {
            for(k=1;k<=D;k++)
                fprintf(ptest,"%f ", data[j][k]);
            fprintf(ptest,"\n");
        }

    }
    fclose(ptest);



    struct samples *p_samples;
    int iter;
    double logprobinit, logprobfinal;

    p_samples = loadSample("a.data");
    ///////
    //
    //
    // remember to define DIM
    //
    //
    ///////
    InitHMM(phmm, 5, 0, 1764, p_samples, 0);
    PrintHMM("test_pre.hmm", phmm);
    BaumWelch(phmm, p_samples, &iter, &logprobinit, &logprobfinal, 200);
    PrintHMM("test.hmm", phmm);
    printf("Test\n");
    return 0;
}

HMM* createHMM(int N, int D)
{
    int i,j,k;
        double tmp;
    HMM* phmm = (HMM*) malloc(sizeof(HMM));
    phmm->N = N;
    phmm->D = D;
    printf("sssss\n");

    phmm->pi = (double *) dvector(1, N);
    for(i=1;i<=N;i++)
            phmm->pi[i] = NAN;
    phmm->pi[1] = 0.0;

    phmm->A = (double **) dmatrix(1,N,1,N);
    for(i=1;i<=N;i++)
    {
        for(j=1;j<=N;j++)
            phmm->A[i][j] = NAN;
    }
    for(i=1;i<N;i++)
    {
        tmp = getRand();
        printf("tmp = %f\n", tmp);
        phmm->A[i][i] = log(tmp);
        phmm->A[i][i+1] = log(1-tmp);
    }
    phmm->A[N][N] = 0.0;
    printf("sssss\n");

    phmm->A_topo = (int **) imatrix(1,N,1,N);
    clear_imatrix(phmm->A_topo,1,N,1,N);
    for(i=1;i<=N;i++)
    {
        for(j=1;j<=N;j++)
        {
            if(j==i || j==i+1)
                    phmm->A_topo[i][j] = 1;
            else
                    phmm->A_topo[i][j] = 0;
        }
    }
    printf("sssss\n");

    phmm->B.miu = (double **) dmatrix(1,N,1,D);
    phmm->B.cov = (double **) dmatrix(1,N,1,D);
    phmm->B.cov_inv = (double **) dmatrix(1,N,1,D);
    for(i=1;i<=N;i++)
    {
        for(j=1;j<=D;j++)
        {
            phmm->B.miu[i][j] = getRand()*2-1;
            phmm->B.cov[i][j] = getRand()*0.06;
            phmm->B.cov_inv[i][j] = 1.0/phmm->B.cov[i][j];
        }
    }
    printf("sssss\n");
    return phmm;
}
