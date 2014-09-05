#include "list.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "sample.h"
#include "nrutil.h"

struct samples * loadSample( const char* filename )
{
    int i,j;
    int sample_count = 0;
    FILE *pf = fopen( filename, "rt" );
    if(pf==NULL)
    {
        printf("Error open file!\n");
        return 0;
    }

    fscanf(pf,"sample_count= %d\n",&sample_count);
    
    struct samples * p_samples = (struct samples*)malloc(sizeof(struct samples));
    p_samples->data = 0;
    p_samples->sample_count = sample_count;
    p_samples->feature_count_per_sample = (int *)calloc( sample_count, sizeof(int) );

    struct sample_node *sample_head = (struct sample_node *)malloc(sizeof(struct sample_node));
    LIST_HEAD_INIT(sample_head->list);

    p_samples->sample_head = sample_head;

    int count = 0;
    p_samples->feature_count_max = 0;
    int sample_index = 0;
    while(sample_index<sample_count)
    {
        fscanf( pf, "# %d\n", &count );
        p_samples->feature_count_per_sample[sample_index++] = count;
        if(count>p_samples->feature_count_max)
            p_samples->feature_count_max = count;
        struct sample_node *sample_tmp = (struct sample_node *)malloc(sizeof(sample_node));
        LIST_HEAD_INIT(sample_tmp->feature_head.list);
        while(count-->0)
        {
            struct feature_node *feature_tmp = (struct feature_node *)malloc(sizeof(struct feature_node));
            for(i=0;i<DIM;i++)
            {
                fscanf(pf, "%lf", &feature_tmp->feature[i]);
            }
            list_add_tail(&feature_tmp->list, &sample_tmp->feature_head.list);
            fscanf(pf,"\n");
        }
        list_add_tail(&sample_tmp->list, &sample_head->list);
    }
    rebuildSamplesData(p_samples);
    return p_samples;
}

void rebuildSamplesData(struct samples* p_samples)
{
    int i = 0;
    assert(p_samples->data == 0);
    p_samples->data = (double ***) calloc(p_samples->sample_count, sizeof(double**)); 
    struct sample_node *sample_head = p_samples->sample_head;

    struct list_head *pos_sample, *pos_feature;
    struct sample_node *p_sample;
    struct feature_node *p_feature;

    int sample_index = 0;
    list_for_each( pos_sample, &sample_head->list )
    {
        p_samples->data[sample_index] = (double**) calloc(p_samples->feature_count_per_sample[sample_index],sizeof(double*));
        p_sample = list_entry( pos_sample, struct sample_node, list );
        int feature_index = 0;
        list_for_each( pos_feature, &p_sample->feature_head.list )
        {
            p_feature = list_entry( pos_feature, struct feature_node, list );
            p_samples->data[sample_index][feature_index] = p_feature->feature;
            feature_index++;
        }
        sample_index++;
    }
    
}

void printSample(struct samples * p_samples)
{
    int i = 0;
    struct sample_node *sample_head = p_samples->sample_head;

    struct list_head *pos_sample, *pos_feature;
    struct sample_node *p_sample;
    struct feature_node *p_feature;
    printf("sample_count= %d\n", p_samples->sample_count);

    int sample_index = 0;
    list_for_each( pos_sample, &sample_head->list )
    {
        printf("# %d\n", p_samples->feature_count_per_sample[sample_index++]);
        p_sample = list_entry( pos_sample, struct sample_node, list );
        list_for_each( pos_feature, &p_sample->feature_head.list )
        {
            p_feature = list_entry( pos_feature, struct feature_node, list );
            for(i=0;i<DIM;i++)
            {
                printf("%f ", p_feature->feature[i]);
            }
            printf("\n");
        }
    }
}

void calcMC(struct samples * p_samples, int N, double **miu, double **cov)
{
    int i, j;
    struct sample_node *sample_head = p_samples->sample_head;

    struct list_head *pos_sample, *pos_feature;
    struct sample_node *p_sample;
    struct feature_node *p_feature;

    int *index_count = (int*)calloc(N,sizeof(int));
    for(i=0;i<N;i++)
        index_count[i] = 0;

    int sample_index = 0;
    list_for_each( pos_sample, &sample_head->list )
    {
        int feature_count = p_samples->feature_count_per_sample[sample_index];
        p_sample = list_entry( pos_sample, struct sample_node, list );
        int feature_index = 0;
        list_for_each( pos_feature, &p_sample->feature_head.list )
        {
            p_feature = list_entry( pos_feature, struct feature_node, list );
            int index = feature_index * N / feature_count + 1;
            index_count[index]++;
            for(i=0;i<DIM;i++)
            {
                miu[index][i+1] = p_feature->feature[i];
                cov[index][i+1] = p_feature->feature[i] * p_feature->feature[i];
            }
            feature_index++;
        }
        sample_index++;
    }

    for(i=1;i<=N;i++)
    {
        for(j=1;j<=DIM;j++)
        {
            miu[i][j] /= index_count[i-1];
            cov[i][j] /= index_count[i-1];
            cov[i][j] -= miu[i][j] * miu[i][j];
        }
    }
    free(index_count);
}

