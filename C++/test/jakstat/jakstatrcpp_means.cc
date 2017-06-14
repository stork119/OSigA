/*
 * -----------------------------------------------------------------
 *  * -----------------------------------------------------------------
 */

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

/* user */
/* #include "model.h" */
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


/* Problem Constants */

#define NEQ   17                /* number of equations  */
#define NPAR  11                /* number of equations  */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(1.0)      /* first output time      */
#define TMULT RCONST(1.0)     /* output time factor     */
#define NOUT  100               /* number of output times */
#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL RCONST(1.0e-4) 

/* Functions Called by the Solver */

struct userData {
    double* parameters; 
    double  stm;
};



struct solverData {
    std::vector< std::vector<double> > trajectory;  
    int solver_flag;
};


static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);


static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* #include "stimulus.h" */

int findNearestNeighbourIndex( double value, double *x, int len );

/*void interp1(double *x, int x_tam, double *y,
  double *xx, int xx_tam, double *yy);
*/
double cosineInterpolation(double x1, double x2,
  double y1, double y2, const realtype t);

double stimulus(const realtype t);



/* Private functions to output5results */
static void writeToFile(N_Vector y, FILE* outputFile);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);



/* pthread */ 

struct threadArgs {
    double* parameters; 
    double* variables;
    std::vector<double> tmesh;
    double stm;
    std::vector< std::vector<double> > output;
    sig_atomic_t *cancel_flag; 
    int solver_flag;
};

/*
 *---------------------------------
 * Run odesolver
 *---------------------------------
 */
solverData run_solver(
  double* parameters, double* variables, double stm)
{
  struct solverData solver_results;
  solver_results.solver_flag = 0;
  std::vector< std::vector<double> > results;
  realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag, iout;

  y = abstol = NULL;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)){
    return(solver_results);
  }
  abstol = N_VNew_Serial(NEQ); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)){
    return(solver_results);
  }
  
  std::vector<double> yvector;  
  for(int i = 0; i < NEQ; ++i){
    double yvar, atol;
    Ith(y,i+1) = variables[i + 1];
    Ith(abstol,i+1) = ATOL;
    yvector.push_back(variables[i + 1]);
  } 
  results.push_back(yvector);

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)){
    return(solver_results);
  }

  /* Call CVODESetuserdata */
  struct userData user_data;
  user_data.parameters = parameters;
  user_data.stm = stm;
  flag = CVodeSetUserData(cvode_mem, (void *) &user_data);
  if (check_flag(&flag, "CVodeSetUserData", 1)){
    return(solver_results);
  }
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)){
    return(solver_results);
  }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)){
    return(solver_results);
  }

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)){
    return(solver_results);
  }

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)){
    return(solver_results);
  }

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  solver_results.solver_flag = 1;
  iout = 0;
  tout = T1;
  while(1) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    std::vector<double> yvector;  
//    printf("%lf\t", t);
    for(int i = 0; i < NEQ; ++i){
      yvector.push_back(Ith(y, i+1));
    //  printf("%lf\t", Ith(y, i+1));
    } 
    results.push_back(yvector);
    //printf("\n");

    if (check_flag(&flag, "CVode", 1)){
      solver_results.solver_flag = 0;
      break;
    }
    if (flag == CV_SUCCESS) {
      iout++;
      tout = tout + TMULT;
    }
    if (iout == NOUT) break;
  }
  //fclose(outputFile);
 
  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free y vector */
  N_VDestroy_Serial(abstol);


   /* Free integrator memory */
  CVodeFree(&cvode_mem);
  solver_results.trajectory = results;
  return(solver_results);
}



void *run_solver_pthread(void* thread_args_ptr){
  struct threadArgs  *thread_args = (struct threadArgs  *) thread_args_ptr;
  struct solverData output = run_solver(
    thread_args->parameters,
    thread_args->variables,
    thread_args->stm);
  thread_args->output = output.trajectory;	
  thread_args->solver_flag = output.solver_flag;	
  *(thread_args->cancel_flag) = 1;	
  pthread_exit((void *) thread_args);

}
 
/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

// [[Rcpp::export]] 
Rcpp::List rmainmean(
  std::vector<double> parameters,
  std::vector<double> variables,
  std::vector<double> tmesh,
  double stm,
  unsigned long time_interval,
  unsigned long time_computation
)
{
//  clock_t start = clock();
  Rcpp::List list;
  list["success"] = 0;

  if( parameters.size() != NPAR ){
    return(list);

  }
  double parametersArray[NPAR + 1];
  for(int i = 0; i < NPAR; ++i){
    parametersArray[i+1] = parameters[i];
  } 

  if( variables.size() < NEQ ){
    return(list);
  }
  double variablesArray[NEQ + 1];
  for(int i = 0; i < NEQ; ++i){
    variablesArray[i+1] = variables[i];
  } 


/*  double tmeshArray[tmesh.size()];
  for(int i = 0; i < tmesh.size(); ++i){
    tmeshArray[i+1] = tmesh[i];
  } 
*/
  /* Run solver */
  sig_atomic_t cancel_flag = 0;

  try
  {	
    struct threadArgs thread_args;
    thread_args.variables = variablesArray;
    thread_args.parameters = parametersArray;
    thread_args.tmesh = tmesh;
    thread_args.stm = stm;
    thread_args.cancel_flag = &cancel_flag;

    pthread_t thread;
    pthread_create(&thread, NULL, run_solver_pthread, (void *) &thread_args);

    void *status = (void*) -1;

    for(int i = 0; i < (int) (time_computation/time_interval); ++i)
    {
      if(cancel_flag != 1){
        usleep(time_interval*1000);
      } else {
        pthread_join(thread, &status);
        break;
      }
    }

    if(cancel_flag != 1){
      pthread_cancel(thread);
      list["success"] = 0;
      list["message"] = "Timeout";
    } else if(thread_args.solver_flag != 1){
      list["success"] = 0;
      list["message"] = "Solver error";
    } else {
	list["success"] = 1;
	list["output"] = thread_args.output;
    }
  } catch (char * message) {
	std::string message_str =  message; 
	list["message"] =   message_str;
        list["success"] = 0;
   }

  return list;

}
/*
int main(){
  std::vector<double> parameters;
  FILE *parFile;
  char *parFilePath = "par.txt";
  parFile = fopen(parFilePath, "r");  
  for(int i = 0; i < NPAR; ++i){
    double tmp;
    if(fscanf(parFile, "%lf", &tmp) != 1){
      return(1);
    }
//    printf("%lf\n",tmp);
    parameters.push_back(tmp);
  } 
  fclose(parFile);

  std::vector<double> variables;
  FILE *varFile;
  char *varFilePath = "var.txt"; 
  varFile = fopen(varFilePath, "r");  
  for(int i = 0; i < NEQ; ++i){
    double yvar;
    if(fscanf(varFile, "%lf", &yvar) != 1){
      return(1);
    }
//    printf("%lf\n",yvar);
    variables.push_back(yvar);
  } 
  fclose(varFile);
   
  rmain(parameters, variables);  

  return(0);

}*/

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */



/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void writeToFile(N_Vector y, FILE* outputFile)
{
  for(int i = 1; i < NEQ; ++i)
  {
    fprintf(outputFile, "%lf,", Ith(y, i) );
  }
  fprintf(outputFile, "%lf\n", Ith(y, NEQ));
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}



/* stimulus.cc */

int findNearestNeighbourIndex( double value, double *x, int len )
{
    double dist;
    int idx;
    int i;

    idx = -1;
    dist = DBL_MAX;
    for ( i = 0; i < len; i++ ) {
        double newDist = value - x[i];
        if ( newDist > 0 && newDist < dist ) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

/*void interp1(double *x, int x_tam, double *y, double *xx, int xx_tam, double *yy)
{
    double dx, dy, *slope, *intercept;
    int i, indiceEnVector;

    slope=(double *)calloc(x_tam,sizeof(double));
    intercept=(double *)calloc(x_tam,sizeof(double));

    for(i = 0; i < x_tam; i++){
        if(i<x_tam-1){
            dx = x[i + 1] - x[i];
            dy = y[i + 1] - y[i];
            slope[i] = dy / dx;
            intercept[i] = y[i] - x[i] * slope[i];
        }else{
            slope[i]=slope[i-1];
            intercept[i]=intercept[i-1];
        }
    }

    for (i = 0; i < xx_tam; i++) {
        indiceEnVector = findNearestNeighbourIndex( xx[i], x, x_tam);
        if (indiceEnVector != -1)
            yy[i]=slope[indiceEnVector] * xx[i] + intercept[indiceEnVector];
        else
            yy[i]=DBL_MAX;
    }
    free(slope);
    free(intercept);
}*/

double cosineInterpolation(double x1, double x2,
  double y1, double y2, const realtype t){
  
  double mu = (t-x1)/(x2-x1);
  double t2 = (1 - cos(mu*M_PI))/2;
  double t1 = y1*(1-t2)+y2*t2;

  return(t1);

}

double stimulus(const realtype t, const realtype stm_val){
  /*if(t > 0 && t < 16){ 
    double yy[] = {0};
    int stm_length = 100;
    double x[stm_length];
    double yi[];
    double xx[]={t};
    int x_tam= sizeof(x)/sizeof(x[0]); 
    interp1(x, x_tam, yi, xx, 1, yy);
    printif("scanf %lf \t %lf \n",t, yy[0]);
    return yy[0];
  } */
  double stm = 0.0;
  if(t >= 0 && t < 0.01){
    stm = cosineInterpolation(0,0.01,0,1,t);
  } else if( t >= 0.01 && t <= 4.99){
     stm = stm_val;  
  } else if( t > 4.99 && t < 5){
     stm = cosineInterpolation(4.99,5,1,0,t);
  }
 // printf("%lf \t %lf \n", t, stm);
  return(stm);
//  return(1);
}

/*void stimulus_test(char *folderPath)
{
  FILE *interOutputFile;
  char interFilePath[100];
  strcpy(interFilePath, folderPath);
  strcat(interFilePath, "/");
  strcat(interFilePath, "interpolation.csv");
  interOutputFile = fopen(interFilePath, "w+");  
  for(realtype t = 0.0; t < 5.0; t += 0.01){
    fprintf(interOutputFile, "%lf,", stimulus(t));
  }                    
  fprintf(interOutputFile, "%lf", stimulus(5.0));
  fclose(interOutputFile);
}*/



/* model.cc */

/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype t, N_Vector yith, N_Vector ydot, void *_user_data)
{
  realtype y[NEQ + 1]; 
  
  struct userData* user_data = (struct userData *) _user_data; 
  double *p = user_data->parameters;
  double stm = stimulus(t, user_data->stm);

  for(int i = 1; i <=  NEQ; ++i){
    y[i] = Ith(yith,i); 
  }

/*
Ith(ydot, 1) =  p[4]*y[13]*2.0-((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 2) =  p[2]*(y[2]*y[2])*-2.0+((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 3) =  -p[3]*y[3]+p[2]*(y[2]*y[2]);
Ith(ydot, 4) =  p[3]*y[3]-p[4]*y[4];
Ith(ydot, 5) =  p[4]*y[4]-p[4]*y[5];
Ith(ydot, 6) =  p[4]*y[5]-p[4]*y[6];
Ith(ydot, 7) =  p[4]*y[6]-p[4]*y[7];
Ith(ydot, 8) =  p[4]*y[7]-p[4]*y[8];
Ith(ydot, 9) =  p[4]*y[8]-p[4]*y[9];
Ith(ydot, 10) =  p[4]*y[9]-p[4]*y[10];
Ith(ydot, 11) =  p[4]*y[10]-p[4]*y[11];
Ith(ydot, 12) =  p[4]*y[11]-p[4]*y[12];
Ith(ydot, 13) =  p[4]*y[12]-p[4]*y[13];
Ith(ydot, 14) =  p[3]*y[3]-p[4]*y[13];
Ith(ydot, 15) =  p[8]*y[16]-stm*p[5]*p[7]*y[15];
Ith(ydot, 16) =  -p[8]*y[16]-p[9]*y[16]+stm*p[5]*p[7]*y[15];
Ith(ydot, 17) =  p[9]*y[16]-p[10]*y[17];
*/

/*hill receptor */
Ith(ydot, 1) =  p[4]*y[13]*2.0-((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 2) =  p[2]*(y[2]*y[2])*-2.0+((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 3) =  -p[3]*y[3]+p[2]*(y[2]*y[2]);
Ith(ydot, 4) =  p[3]*y[3]-p[4]*y[4];
Ith(ydot, 5) =  p[4]*y[4]-p[4]*y[5];
Ith(ydot, 6) =  p[4]*y[5]-p[4]*y[6];
Ith(ydot, 7) =  p[4]*y[6]-p[4]*y[7];
Ith(ydot, 8) =  p[4]*y[7]-p[4]*y[8];
Ith(ydot, 9) =  p[4]*y[8]-p[4]*y[9];
Ith(ydot, 10) =  p[4]*y[9]-p[4]*y[10];
Ith(ydot, 11) =  p[4]*y[10]-p[4]*y[11];
Ith(ydot, 12) =  p[4]*y[11]-p[4]*y[12];
Ith(ydot, 13) =  p[4]*y[12]-p[4]*y[13];
Ith(ydot, 14) =  p[3]*y[3]-p[4]*y[13];
Ith(ydot, 15) =  p[8]*y[16]-(stm*p[5]*p[7]*(y[15]*y[15]))/(y[15]*y[15]+p[11]);
Ith(ydot, 16) =  -p[8]*y[16]-p[9]*y[16]+(stm*p[5]*p[7]*(y[15]*y[15]))/(y[15]*y[15]+p[11]);
Ith(ydot, 17) =  p[9]*y[16]-p[10]*y[17];

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(long int N, realtype t,
               N_Vector yvec, N_Vector fy, DlsMat J, void *_user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y[NEQ + 1];

  struct userData* user_data = (struct userData *) _user_data; 
  double *p = user_data->parameters;
  double stm = stimulus(t, user_data->stm);

  for(int i = 1; i <= NEQ; ++i){
    y[i] = Ith(yvec,i); 
  }

  for(int i = 1; i <= NEQ; ++i){
    for(int j = 1; j <= NEQ; ++j){
       IJth(J,i,j) = RCONST(0.0);
    }
  }
  
/*
IJth(J, 1, 1) =  -((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 1, 13) =  p[4]*2.0;
IJth(J, 1, 16) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 1, 17) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 2, 1) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 2, 2) =  p[2]*y[2]*-4.0;
IJth(J, 2, 16) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 2, 17) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 3, 2) =  p[2]*y[2]*2.0;
IJth(J, 3, 3) =  -p[3];
IJth(J, 4, 3) =  p[3];
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
IJth(J, 14, 3) =  p[3];
IJth(J, 14, 13) =  -p[4];
IJth(J, 15, 15) =  -stm*p[5]*p[7];
IJth(J, 15, 16) =  p[8];
IJth(J, 16, 15) =  stm*p[5]*p[7];
IJth(J, 16, 16) =  -p[8]-p[9];
IJth(J, 17, 16) =  p[9];
IJth(J, 17, 17) =  -p[10];
*/

/*hill receptor */
IJth(J, 1, 1) =  -((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 1, 13) =  p[4]*2.0;
IJth(J, 1, 16) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 1, 17) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 2, 1) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 2, 2) =  p[2]*y[2]*-4.0;
IJth(J, 2, 16) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 2, 17) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 3, 2) =  p[2]*y[2]*2.0;
IJth(J, 3, 3) =  -p[3];
IJth(J, 4, 3) =  p[3];
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
IJth(J, 14, 3) =  p[3];
IJth(J, 14, 13) =  -p[4];
IJth(J, 15, 15) =  (stm*p[5]*p[7]*y[15]*-2.0)/(y[15]*y[15]+p[11])+stm*1.0/pow(y[15]*y[15]+p[11],2.0)*p[5]*p[7]*(y[15]*y[15]*y[15])*2.0;
IJth(J, 15, 16) =  p[8];
IJth(J, 16, 15) =  (stm*p[5]*p[7]*y[15]*2.0)/(y[15]*y[15]+p[11])-stm*1.0/pow(y[15]*y[15]+p[11],2.0)*p[5]*p[7]*(y[15]*y[15]*y[15])*2.0;
IJth(J, 16, 16) =  -p[8]-p[9];
IJth(J, 17, 16) =  p[9];
IJth(J, 17, 17) =  -p[10];


  return(0);
}


