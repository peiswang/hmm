/*
**      File:   nrutil.h
**      Purpose: Memory allocation routines borrowed from the
**              book "Numerical Recipes" by Press, Flannery, Teukolsky,
**              and Vetterling.
**              state sequence and probablity of observing a sequence
**              given the model.
**      Organization: University of Maryland
**
**      $Id: nrutil.h,v 1.2 1998/02/19 16:32:42 kanungo Exp kanungo $
*/

#ifndef __NRUTIL_H__
#define __NRUTIL_H__

#include <math.h>


double *dvector(int,int);
int *ivector(int,int);
double **dmatrix(int,int,int,int);
int **imatrix(int,int,int,int);
double ***dmatrix3d(int,int,int,int,int,int);
double ****dmatrix4d(int,int,int,int,int,int,int,int);

void clear_dvector(double*,int,int);
void clear_dvector_nan(double*,int,int);
void clear_ivector(int*,int,int);
void clear_dmatrix(double**,int,int,int,int);
void clear_dmatrix_nan(double**,int,int,int,int);
void clear_imatrix(int**,int,int,int,int);
void clear_dmatrix3d(double***,int,int,int,int,int,int);

void free_dvector(double*,int,int);
void free_ivector(int*,int,int);
void free_dmatrix(double**,int,int,int,int);
void free_imatrix(int**,int,int,int,int);
void free_dmatrix3d(double***,int,int,int,int,int,int);
void free_dmatrix4d(double****,int,int,int,int,int,int,int,int);;

double logSum(double x, double y);
double logExp(double x);


//double logProd(double x, double y)
//{
//    return x + y;
//}
//
//inline double logProd3(double x, double y, double z)
//{
//    if(isnan(x) || isnan(y) || isnan(z))
//        return NAN;
//    return x + y + z;
//}
//
//inline double logProd4(double x, double y, double z, double w)
//{
//    if(isnan(x) || isnan(y) || isnan(z) || isnan(w))
//        return NAN;
//    return x + y + z + w;
//}
//
//inline double logProd5(double x, double y, double z, double w, double v)
//{
//    if(isnan(x) || isnan(y) || isnan(z) || isnan(w) || isnan(v))
//        return NAN;
//    return x + y + z + w + v;
//}


#endif
