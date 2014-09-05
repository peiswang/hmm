#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include "list.h"
#include <stdio.h>
#include <stdlib.h>

#define DIM 1764 

typedef struct feature_node
{
    double feature[DIM];
    struct list_head list;
} feature_node;

typedef struct sample_node
{
    struct feature_node feature_head;
    struct list_head list;
} sample_node;

struct samples
{
    struct sample_node * sample_head;
    int sample_count;
    int feature_count_max;
    int * feature_count_per_sample;
    double ***data;
};

struct samples * loadSample( const char* filename );

void rebuildSamplesData(struct samples *);

void printSample(struct samples * p_samples);

#endif
