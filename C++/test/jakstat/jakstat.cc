/*
 * -----------------------------------------------------------------
 * $Revision: 4834 $
 * $Date: 2016-08-01 16:59:05 -0700 (Mon, 01 Aug 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h> 
#include <math.h>
#include <string.h>

/* Header files with a description of contents used */

#include <cvodes/cvodes.h>           /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

/* User-defined vector and matrix accessor macros: Ith, IJth */

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

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   17                /* number of equations  */
#define NPAR  10                /* number of equations  */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(5.0)      /* first output time      */
#define TMULT RCONST(5.0)     /* output time factor     */
#define NOUT  20               /* number of output times */
#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL RCONST(1.0e-4) 

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);


static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void writeToFile(N_Vector y, FILE* outputFile);
static void PrintOutput(realtype t, realtype y1);
static void PrintRootInfo(int root_f1, int root_f2);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);


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

void interp1(double *x, int x_tam, double *y, double *xx, int xx_tam, double *yy)
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
}

double cosineInterpolation(double x1, double x2, double y1, double y2, const realtype t){
  
  double mu = (t-x1)/(x2-x1);
  double t2 = (1 - cos(mu*M_PI))/2;
  double t1 = y1*(1-t2)+y2*t2;

  return(t1);

}

double stimulus(const realtype t){
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
     stm = 1;  
  } else if( t > 4.99 && t < 5){
     stm = cosineInterpolation(4.99,5,1,0,t);
  }
 // printf("%lf \t %lf \n", t, stm);
  return(stm);
//  return(1);
}


void stimulus_test(char *folderPath)
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
}

/*
 *---------------------------------
 * Run odesolver
 *---------------------------------
 */
int run_solver(char *folderPath)
{
  realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag, flagr, iout;
  int rootsfound[2];

  y = abstol = NULL;
  cvode_mem = NULL;

  /* Read parameters from file */ 
  double parameters[NPAR + 1];
  FILE *parFile;
  char *parFilePath = "par.txt";
  parFile = fopen(parFilePath, "r");  
  for(int i = 0; i < NPAR; ++i){
    if(fscanf(parFile, "%lf", &parameters[i+1]) != 1){
      return(1);
    }
  } 
  fclose(parFile);
  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(NEQ); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);
  
  /* Read initial variables form file var.txt */
  FILE *varFile;
  char *varFilePath = "var.txt"; 
  varFile = fopen(varFilePath, "r");  
  for(int i = 0; i < NEQ; ++i){
    double yvar, atol;
    if(fscanf(varFile, "%lf %lf", &yvar, &atol) != 2){
      return(1);
    }
    Ith(y,i+1) = yvar;
    Ith(abstol,i+1) = ATOL;
  } 
  fclose(varFile);

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Initialize output file */
  FILE *outputFile;
  char outFilePath[100];
  strcpy(outFilePath, folderPath);
  strcat(outFilePath, "/");
  strcat(outFilePath, "out.csv");
  outputFile = fopen(outFilePath, "w+");  
  
  writeToFile(y,outputFile);
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  
  /* Call CVODESetuserdata */
  flag = CVodeSetUserData(cvode_mem, (void *) parameters);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  /* Call CVodeRootInit to specify the root function g with 2 components */
/*  flag = CVodeRootInit(cvode_mem, 2, g);
  if (check_flag(&flag, "CVodeRootInit", 1)) return(1);*/

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  //printf(" \n3-species kinetics problem\n\n");
  iout = 0;
  tout = T1;
  t = T0;
  while(1) {
    //printf("KAROL \t %lf\n", tout);
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    //PrintOutput(t, Ith(y,13));
    writeToFile(y,outputFile);
    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      iout++;
      tout = tout + TMULT;
    }

    if (iout == NOUT) break;
  }
  fclose(outputFile);
  /* Print some final statistics */
  //  PrintFinalStats(cvode_mem);
 
  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free y vector */
  N_VDestroy_Serial(abstol);


   /* Free integrator memory */
  CVodeFree(&cvode_mem);
 

  return(0);
}
 
/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, const char *argv[])
{
  clock_t start = clock();
  
  for(int run_i = 0; run_i < 1000; ++run_i)
  { 
    int flag;
    
    /* Creat output folder */
    char folderPath[100];
    strcpy(folderPath, "test1");
    flag = mkdir(folderPath, 0755);
    if (flag == -1) 
    {
//      perror(argv[0]);
      // exit(EXIT_FAILURE);
    }
    strcat(folderPath, "/stm-normal-"); 
    char run_str_i[10];
    sprintf(run_str_i, "%d", run_i);
    strcat(folderPath,  run_str_i);  
    //printf("%s\n", folderPath);  
    flag = mkdir(folderPath, 0755);
    if (flag == -1) 
    {
//      perror(argv[0]);
      // exit(EXIT_FAILURE);
    }

    /* Run solver */
    if(run_solver(folderPath) != 0){ return(1); }
    
  }
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;
  printf("%f\n", seconds);
  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */



/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype t, N_Vector yith, N_Vector ydot, void *user_data)
{
  realtype y[NEQ + 1]; 

  double *p = (double *) user_data;
  double stm = stimulus(t);

  for(int i = 1; i <=  NEQ; ++i){
    y[i] = Ith(yith,i); 
  }
  //printf("y17\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, stm, Ith(ydot,17),  p[9], y[16], p[10], y[17]);  

  Ith(ydot,1) = 2*p[4]*y[13] - (p[1]*y[1]*(y[16] + y[17]))/(p[6] + y[1]);
  Ith(ydot,2) = (p[1]*y[1]*(y[16] + y[17]))/(p[6] + y[1]) - 2*p[2]*y[2]*y[2];
  Ith(ydot,3) = p[2]*y[2]*y[2] - p[3]*y[3];
  Ith(ydot,4) = p[3]*y[3] - p[4]*y[4];
  Ith(ydot,5) = p[4]*y[4] - p[4]*y[5];
  Ith(ydot,6) = p[4]*y[5] - p[4]*y[6]; 
  Ith(ydot,7) = p[4]*y[6] - p[4]*y[7]; 
  Ith(ydot,8) = p[4]*y[7] - p[4]*y[8]; 
  Ith(ydot,9) = p[4]*y[8] - p[4]*y[9]; 
  Ith(ydot,10) = p[4]*y[9] - p[4]*y[10];
  Ith(ydot,11) = p[4]*y[10] - p[4]*y[11];
  Ith(ydot,12) = p[4]*y[11] - p[4]*y[12];
  Ith(ydot,13) = p[4]*y[12] - p[4]*y[13];
  Ith(ydot,14) = p[3]*y[3] - p[4]*y[13];
  Ith(ydot,15) = p[8]*y[16] - p[5]*p[7]*y[15]*stm;
  Ith(ydot,16) = p[5]*p[7]*y[15]*stm - p[9]*y[16] - p[8]*y[16];
  Ith(ydot,17) = p[9]*y[16] - p[10]*y[17];

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(long int N, realtype t,
               N_Vector yvec, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y[NEQ + 1];

  double *p = (double *) user_data;
  double stm = stimulus(t);

  for(int i = 1; i <= NEQ; ++i){
    y[i] = Ith(yvec,i); 
  }

  for(int i = 1; i <= NEQ; ++i){
    for(int j = 1; j <= NEQ; ++j){
       IJth(J,i,j) = RCONST(0.0);
    }
  }
  
  IJth(J, 1, 1) = (p[1]*y[1]*(y[16] + y[17]))/((p[6] + y[1])*(p[6] + y[1]))  - (p[1]*(y[16] + y[17]))/(p[6] + y[1]);
  IJth(J, 1, 13) = 2*p[4]; 
  IJth(J, 1, 17) =  -(p[1]*y[1])/(p[6] + y[1]);
  IJth(J, 1, 16) = -(p[1]*y[1])/(p[6] + y[1]);
  IJth(J, 2, 1) = (p[1]*(y[16] + y[17]))/(p[6] + y[1]) - (p[1]*y[1]*(y[16] + y[17]))/((p[6] + y[1])*(p[6] + y[1]));
  IJth(J, 2, 2) = -4*p[2]*y[2];
  IJth(J, 2, 16) = p[1]*y[1]/(p[6] + y[1]);
  IJth(J, 2, 17) = p[1]*y[1]/(p[6] + y[1]);
  IJth(J, 3, 2) = 2*p[2]*y[2]; 
  IJth(J, 3, 3) = -p[3];
  IJth(J, 4, 3) = p[3];
  IJth(J, 4, 4) = -p[4];
  IJth(J, 5, 4) = p[4];
  IJth(J, 5, 5) = -p[4];
  IJth(J, 6, 5) = p[4];
  IJth(J, 6, 6) = -p[4];
  IJth(J, 7, 6) = p[4];
  IJth(J, 7, 7) = -p[4];
  IJth(J, 8, 7) = p[4];
  IJth(J, 8, 8) = -p[4];
  IJth(J, 9, 8) = p[4];
  IJth(J, 9, 9) = -p[4];
  IJth(J, 10, 9) = p[4];
  IJth(J, 10, 10) = -p[4];
  IJth(J, 11, 10) = p[4];
  IJth(J, 11, 11) = -p[4];
  IJth(J, 12, 11) = p[4];
  IJth(J, 12, 12) = -p[4];
  IJth(J, 13, 12) = p[4];
  IJth(J, 13, 13) = -p[4];
  IJth(J, 14, 3) = p[3];
  IJth(J, 14, 13) = -p[4];
  IJth(J, 15, 15) = -p[5]*p[7]*stm;
  IJth(J, 15, 16) = p[8];
  IJth(J, 16, 15) = p[5]*p[7]*stm;
  IJth(J, 16, 16) = - p[8] - p[9];
  IJth(J, 17, 16) = p[9];
  IJth(J, 17, 17) = -p[10];
 
  return(0);
}

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

static void PrintOutput(realtype t, realtype y1)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le\n", t, y1);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %lf      y =%lf\n", t, y1);
#else
  printf("At t = %0.4e      y =%14.6e\n", t, y1);
#endif

  return;
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
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
