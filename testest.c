#include <stdio.h>
#include "hmm.h"
#include "sample.h"

int main()
{
    HMM hmm;
    HMM *phmm = &hmm;
    struct samples *p_samples;
    int iter;
    double logprobinit, logprobfinal;

    p_samples = loadSample("test.data");
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
