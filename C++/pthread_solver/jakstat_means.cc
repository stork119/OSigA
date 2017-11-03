#include "jakstat_means.h"

#define NEQ_jakstat_means   8                /* number of equations  */
#define NPAR_jakstat_means  14                /* number of equations  */

struct modelData returnModelStructure_jakstat_means(){
  struct modelData model_data;
  model_data.neq = NEQ_jakstat_means;
  model_data.npar = NPAR_jakstat_means;
  return(model_data); 
}

/*
 * f routine. Compute function f(t,y). 
 */

int f_jakstat_means(
  realtype t,
  N_Vector yith, 
  N_Vector ydot,
  void *_user_data)
{
  realtype y[NEQ_jakstat_means + 1]; 
  
  struct userData* user_data = (struct userData *) _user_data; 
  double *p = user_data->parameters;
  double stm = stimulus(t, user_data->stm);

  for(int i = 1; i <=  NEQ_jakstat_means; ++i){
    y[i] = Ith(yith,i); 
  }

  Ith(ydot, 1) =  (p[12]*p[14]*y[5])/p[11]-((p[11]*y[8]+p[13]*y[7])*p[1]*y[1])/((p[6]+y[1])*p[11]);
  Ith(ydot, 2) =  p[2]*(y[2]*y[2])*-2.0+((p[11]*y[8]+p[13]*y[7])*p[1]*y[1])/((p[6]+y[1])*p[11]); 
  Ith(ydot, 3) =  -p[3]*y[3]+p[2]*(y[2]*y[2]);
  Ith(ydot, 4) =  -p[4]*y[4]+(p[3]*p[11]*y[3])/p[12];
  Ith(ydot, 5) =  p[4]*y[4]*2.0-p[14]*y[5];
  Ith(ydot, 6) =  p[8]*y[7]-stm*p[5]*p[7]*y[6];
  Ith(ydot, 7) =  -p[8]*y[7]-p[9]*y[7]+stm*p[5]*p[7]*y[6];
  Ith(ydot, 8) =  -p[10]*y[8]+(p[9]*p[13]*y[7])/p[11];

  /*Ith(ydot, 1) =  (p[4]*p[12]*y[13]*2.0)/p[11]-((p[11]*y[17]+p[13]*y[16])*p[1]*y[1])/(p[11]*y[1]+p[6]);
  Ith(ydot, 2) =  p[2]*p[11]*(y[2]*y[2])*-2.0+((p[11]*y[17]+p[13]*y[16])*p[1]*y[1])/(p[11]*y[1]+p[6]);
  Ith(ydot, 3) =  -p[3]*y[3]+p[2]*p[11]*(y[2]*y[2]);
  Ith(ydot, 4) =  -p[4]*y[4]+(p[3]*p[11]*y[3])/p[12];
  Ith(ydot, 5) =  p[4]*y[4]-p[4]*y[5];
  Ith(ydot, 6) =  p[4]*y[5]-p[4]*y[6];
  Ith(ydot, 7) =  p[4]*y[6]-p[4]*y[7];
  Ith(ydot, 8) =  p[4]*y[7]-p[4]*y[8];
  Ith(ydot, 9) =  p[4]*y[8]-p[4]*y[9];
  Ith(ydot, 10) =  p[4]*y[9]-p[4]*y[10];
  Ith(ydot, 11) =  p[4]*y[10]-p[4]*y[11];
  Ith(ydot, 12) =  p[4]*y[11]-p[4]*y[12];
  Ith(ydot, 13) =  p[4]*y[12]-p[4]*y[13];
  Ith(ydot, 14) =  -p[4]*y[13]+(p[3]*p[11]*y[3])/p[12];
  Ith(ydot, 15) =  p[8]*y[16]-stm*p[5]*p[7]*y[15];
  Ith(ydot, 16) =  -p[8]*y[16]-p[9]*y[16]+stm*p[5]*p[7]*y[15];
  Ith(ydot, 17) =  (p[9]*p[13]*y[16])/p[12]-(p[10]*p[11]*y[17])/p[12];*/

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

int Jac_jakstat_means(long int N, realtype t,
               N_Vector yvec, N_Vector fy, DlsMat J, void *_user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y[NEQ_jakstat_means + 1];

  struct userData* user_data = (struct userData *) _user_data; 
  double *p = user_data->parameters;
  double stm = stimulus(t, user_data->stm);

  for(int i = 1; i <= NEQ_jakstat_means; ++i){ y[i] = Ith(yvec,i); 
  }

  for(int i = 1; i <= NEQ_jakstat_means; ++i){
    for(int j = 1; j <= NEQ_jakstat_means; ++j){
       IJth(J,i,j) = RCONST(0.0);
    }
  }

  IJth(J, 1, 1) =  -((p[11]*y[8]+p[13]*y[7])*p[1])/((p[6]+y[1])*p[11])+(1.0/pow(p[6]+y[1],2.0)*(p[11]*y[8]+p[13]*y[7])*p[1]*y[1])/p[11];
  IJth(J, 1, 5) =  (p[12]*p[14])/p[11];
  IJth(J, 1, 7) =  -(p[1]*p[13]*y[1])/((p[6]+y[1])*p[11]);
  IJth(J, 1, 8) =  -(p[1]*y[1])/(p[6]+y[1]);
  IJth(J, 2, 1) =  ((p[11]*y[8]+p[13]*y[7])*p[1])/((p[6]+y[1])*p[11])-(1.0/pow(p[6]+y[1],2.0)*(p[11]*y[8]+p[13]*y[7])*p[1]*y[1])/p[11];
  IJth(J, 2, 2) =  p[2]*y[2]*-4.0;
  IJth(J, 2, 7) =  (p[1]*p[13]*y[1])/((p[6]+y[1])*p[11]);
  IJth(J, 2, 8) =  (p[1]*y[1])/(p[6]+y[1]);
  IJth(J, 3, 2) =  p[2]*y[2]*2.0;
  IJth(J, 3, 3) =  -p[3];
  IJth(J, 4, 3) =  (p[3]*p[11])/p[12];
  IJth(J, 4, 4) =  -p[4];
  IJth(J, 5, 4) =  p[4]*2.0;
  IJth(J, 5, 5) =  -p[14];
  IJth(J, 6, 6) =  -stm*p[5]*p[7];
  IJth(J, 6, 7) =  p[8];
  IJth(J, 7, 6) =  stm*p[5]*p[7];
  IJth(J, 7, 7) =  -p[8]-p[9];
  IJth(J, 8, 7) =  (p[9]*p[13])/p[11];
  IJth(J, 8, 8) =  -p[10];

/*  IJth(J, 1, 1) =  -((p[11]*y[17]+p[13]*y[16])*p[1])/(p[11]*y[1]+p[6])+1.0/pow(p[11]*y[1]+p[6],2.0)*(p[11]*y[17]+p[13]*y[16])*p[1]*p[11]*y[1];
  IJth(J, 1, 13) =  (p[4]*p[12]*2.0)/p[11];
  IJth(J, 1, 16) =  -(p[1]*p[13]*y[1])/(p[11]*y[1]+p[6]);
  IJth(J, 1, 17) =  -(p[1]*p[11]*y[1])/(p[11]*y[1]+p[6]);
  IJth(J, 2, 1) =  ((p[11]*y[17]+p[13]*y[16])*p[1])/(p[11]*y[1]+p[6])-1.0/pow(p[11]*y[1]+p[6],2.0)*(p[11]*y[17]+p[13]*y[16])*p[1]*p[11]*y[1];
  IJth(J, 2, 2) =  p[2]*p[11]*y[2]*-4.0;
  IJth(J, 2, 16) =  (p[1]*p[13]*y[1])/(p[11]*y[1]+p[6]);
  IJth(J, 2, 17) =  (p[1]*p[11]*y[1])/(p[11]*y[1]+p[6]);
  IJth(J, 3, 2) =  p[2]*p[11]*y[2]*2.0;
  IJth(J, 3, 3) =  -p[3];
  IJth(J, 4, 3) =  (p[3]*p[11])/p[12];
  IJth(J, 4, 4) =  -p[4];
  IJth(J, 5, 4) =  p[4];
  IJth(J, 5, 5) =  -p[4];
  IJth(J, 6, 5) =  p[4];
  IJth(J, 6, 6) =  -p[4];
  IJth(J, 7, 6) =  p[4];
  IJth(J, 7, 7) =  -p[4];
  IJth(J, 8, 7) =  p[4];
  IJth(J, 8, 8) =  -p[4];
  IJth(J, 9, 8) =  p[4];
  IJth(J, 9, 9) =  -p[4];
  IJth(J, 10, 9) =  p[4];
  IJth(J, 10, 10) =  -p[4];
  IJth(J, 11, 10) =  p[4];
  IJth(J, 11, 11) =  -p[4];
  IJth(J, 12, 11) =  p[4];
  IJth(J, 12, 12) =  -p[4];
  IJth(J, 13, 12) =  p[4];
  IJth(J, 13, 13) =  -p[4];
  IJth(J, 14, 3) =  (p[3]*p[11])/p[12];
  IJth(J, 14, 13) =  -p[4];
  IJth(J, 15, 15) =  -stm*p[5]*p[7];
  IJth(J, 15, 16) =  p[8];
  IJth(J, 16, 15) =  stm*p[5]*p[7];
  IJth(J, 16, 16) =  -p[8]-p[9];
  IJth(J, 17, 16) =  (p[9]*p[13])/p[12];
  IJth(J, 17, 17) =  -(p[10]*p[11])/p[12];*/

  return(0);
}
