
#ifndef __RAND_H__
#define __RAND_H__

#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>

int getSeed(void)
{
    return (int) getpid();
}

void setSeed(int seed)
{
    srand(seed);
}

double getRand(void)
{
    return (double)rand()/RAND_MAX;
}

double getNormalRand(double mean, double sigma)
{
    int i;
    double x = 0.0;
    for(i=0;i<12;i++)
    {
        x += (double)rand()/RAND_MAX;
    }
    x -= 6.0;
    x = x * sigma + mean;
    return x;
}

#endif
