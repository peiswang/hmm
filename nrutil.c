/*
**      File:   nrutil.c
**      Purpose: Memory allocation routines borrowed from the
**		book "Numerical Recipes" by Press, Flannery, Teukolsky,
**		and Vetterling. 
**              state sequence and probablity of observing a sequence
**              given the model.
**      Organization: University of Maryland
**
**      $Id: nrutil.c,v 1.2 1998/02/19 16:31:35 kanungo Exp kanungo $
*/

#include <malloc.h>
#include <math.h>
#include <stdio.h>

#include "nrutil.h"

static char rcsid[] = "$Id: nrutil.c,v 1.2 1998/02/19 16:31:35 kanungo Exp kanungo $";


int *ivector(int nl,int nh)
{
	int *v;

	v=(int*)calloc((unsigned) (nh-nl+1),sizeof(int));
	if (!v) printf("allocation failure in ivector()");
	return v-nl;
}

void clear_ivector(int *m,int nl,int nh)
{
    int i;
    for(i=nl;i<=nh;i++)
        m[i] = 0;
}

double *dvector(int nl,int nh)
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1),sizeof(double));
	if (!v) printf("allocation failure in dvector()");
	return v-nl;
}

void clear_dvector(double *m,int nl,int nh)
{
    int i;
    for(i=nl;i<=nh;i++)
        m[i] = 0.0;
}

void clear_dvector_nan(double *m,int nl,int nh)
{
    int i;
    for(i=nl;i<=nh;i++)
        m[i] = NAN;
}

int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	int **m;

	m=(int**) calloc((unsigned) (nrh-nrl+1),sizeof(int*));
	if (!m) printf("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int*) calloc((unsigned) (nch-ncl+1),sizeof(int));
		if (!m[i]) printf("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}


double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	if (!m) printf("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
		if (!m[i]) printf("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}


void clear_imatrix(int **m,int  nrl,int nrh,int ncl,int nch)
{
    int i,j;
    for(i=nrl;i<=nrh;i++)
    {
        for(j=ncl;j<=nch;j++)
            m[i][j] = 0;
    }
}

void clear_dmatrix(double **m,int  nrl,int nrh,int ncl,int nch)
{
    int i,j;
    for(i=nrl;i<=nrh;i++)
    {
        for(j=ncl;j<=nch;j++)
            m[i][j] = 0.0;
    }
}

void clear_dmatrix_nan(double **m,int  nrl,int nrh,int ncl,int nch)
{
    int i,j;
    for(i=nrl;i<=nrh;i++)
    {
        for(j=ncl;j<=nch;j++)
            m[i][j] = NAN;
    }
}

double ***dmatrix3d(int nsl,int nsh,int nrl,int nrh,int ncl,int nch)
{
	int i,j;
	double ***m;

	m=(double ***) calloc((unsigned) (nsh-nsl+1),sizeof(double**));
	if (!m) printf("allocation failure 1 in dmatrix3d()");
	m -= nsl;

    for(j=nsl;j<=nsh;j++)
    {
        m[j] = (double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	    if (!m[j]) printf("allocation failure 2 in dmatrix3d()");
        m[j] -= nrl;
        
	    for(i=nrl;i<=nrh;i++) {
	    	m[j][i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
	    	if (!m[j][i]) printf("allocation failure 3 in dmatrix3d()");
	    	m[j][i] -= ncl;
	    }
    }
	return m;
}

double ****dmatrix4d(int nbl,int nbh,int nsl,int nsh,int nrl,int nrh,int ncl,int nch)
{
	int i,j,k;
	double ****m;

	m=(double ****) calloc((unsigned) (nbh-nbl+1),sizeof(double***));
	if (!m) printf("allocation failure 1 in dmatrix4d()");
	m -= nbl;

    for(k=nbl;k<=nbh;k++)
    {
        m[k] = (double ***) calloc((unsigned) (nsh-nsl+1),sizeof(double**));
	    if (!m[k]) printf("allocation failure 2 in dmatrix4d()");
        m[k] -= nsl;
        for(j=nsl;j<=nsh;j++)
        {
            m[k][j] = (double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	        if (!m[k][j]) printf("allocation failure 3 in dmatrix4d()");
            m[k][j] -= nrl;
            
	        for(i=nrl;i<=nrh;i++) {
	        	m[k][j][i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
	        	if (!m[k][j][i]) printf("allocation failure 4 in dmatrix4d()");
	        	m[k][j][i] -= ncl;
	        }
        }
    }
	return m;
}
void clear_dmatrix3d(double ***m,int nsl,int nsh,int  nrl,int nrh,int ncl,int nch)
{
    int i,j,k;
    for(i=nsl;i<=nsh;i++)
    {
        for(j=nrl;j<=nrh;j++)
            for(k=ncl;k<=nch;k++)
                m[i][j][k] = 0.0;
    }
}

void free_ivector(int *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_dvector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix3d(double ***m,int nsl,int nsh,int nrl,int nrh,int ncl,int nch)
{
	int i,j;
    for(j=nsh;j>=nsl;j--)
    {
        for(i=nrh;i>=nrl;i--)
        {
            free((char*) (m[j][i]+ncl));
        }
        free((char*) (m[j]+nrl));
    }
    free((char*)(m+nsl));
}

void free_dmatrix4d(double ****m,int nbl,int nbh,int nsl,int nsh,int nrl,int nrh,int ncl,int nch)
{
	int i,j,k;
    for(k=nbh;k>=nbl;k--)
    {
        for(j=nsh;j>=nsl;j--)
        {
            for(i=nrh;i>=nrl;i--)
            {
                free((char*) (m[k][j][i]+ncl));
            }
            free((char*) (m[k][j]+nrl));
        }
        free((char*) (m[k]+nsl));
    }
    free((char*)(m+nbl));
}

double logSum(double x, double y)
{
    if(isnan(x))
            return y;
    if(isnan(y))
            return x;
    if(x>y)
        return x + log(1+exp(y-x));
    return y + log(1+exp(x-y));
    //return log(exp(x)+exp(y));
}

double logExp(double x)
{
    if(isnan(x))
            return 0;
    return exp(x);
}
