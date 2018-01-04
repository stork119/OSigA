#ifndef MODEL_RCPP_H
#define MODEL_RCPP_H

/* C */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h> 
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>

/* CPP */

#include <vector>

/* Rcpp */
#include <Rcpp.h>

/* Header files with a description of contents used */

#include <cvodes/cvodes.h>           /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */


/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

/* User-defined vector and matrix accessor macros: Ith, IJth */
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

struct userData {
    double* parameters; 
    double  stm;
};

/* */

struct modelData {
    int neq; 
    int npar;
};


/* #include "stimulus.h" */

int findNearestNeighbourIndex( double value, double *x, int len );

void interp1(double *x, int x_tam, double *y,
  double *xx, int xx_tam, double *yy);

double cosineInterpolation(double x1, double x2,
  double y1, double y2, const realtype t);


double stimulus(const realtype t, const realtype stm_val);

#endif
