#ifndef JAKSTAT_MEANS_H
#define JAKSTAT_MEANS_H

#include "model_rcpp_libraries.h"
struct modelData returnModelStructure_jakstat_means();

int f_jakstat_means(realtype t, N_Vector yith, N_Vector ydot, void *_user_data);

int Jac_jakstat_means(long int N, realtype t,
               N_Vector yvec, N_Vector fy, DlsMat J, void *_user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);



#endif
