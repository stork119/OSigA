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

/* C */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h> 
#include <math.h>
#include <string.h>

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

#define NEQ   170                /* number of equations  */
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

double cosineInterpolation(double x1, double x2,
  double y1, double y2, const realtype t){
  
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

/*
 *---------------------------------
 * Run odesolver
 *---------------------------------
 */
std::vector< std::vector<double> > run_solver(
  double* parameters, double* variables)
{
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
    return(std::vector< std::vector<double> >());
  }
  abstol = N_VNew_Serial(NEQ); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)){
    return(std::vector< std::vector<double> >());
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
    return(std::vector< std::vector<double> >());
  }

  /* Call CVODESetuserdata */
  flag = CVodeSetUserData(cvode_mem, (void *) parameters);
  if (check_flag(&flag, "CVodeSetUserData", 1)){
    return(std::vector< std::vector<double> >());
  }
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)){
    return(std::vector< std::vector<double> >());
  }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)){
    return(std::vector< std::vector<double> >());
  }

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)){
    return(std::vector< std::vector<double> >());
  }

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)){
    return(std::vector< std::vector<double> >());
  }

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  iout = 0;
  tout = T1;
  while(1) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    std::vector<double> yvector;  
    for(int i = 0; i < NEQ; ++i){
      yvector.push_back(Ith(y, i+1));
//      printf("%lf\t", Ith(y, i+1));
    } 
    results.push_back(yvector);
//    printf("\n");

    if (check_flag(&flag, "CVode", 1)) break;
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
 
  return(results);
}
 
/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

// [[Rcpp::export]] 
std::vector< std::vector<double> > rmain(std::vector<double> parameters, std::vector<double> variables)
{
//  clock_t start = clock();
  
  if( parameters.size() != NPAR ){
    return(std::vector< std::vector<double> >());

  }
  double parametersArray[NPAR + 1];
  for(int i = 0; i < NPAR; ++i){
    parametersArray[i+1] = parameters[i];
  } 

  if( variables.size() != NEQ ){
    return(std::vector< std::vector<double> >());
  }
  double variablesArray[NEQ + 1];
  for(int i = 0; i < NEQ; ++i){
    variablesArray[i+1] = variables[i];
  } 

    
  /* Run solver */
  std::vector< std::vector<double> > result = run_solver(parametersArray, variablesArray);
    
  /*clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;
  printf("%f\n", seconds);*/
  return(result);
}

/*int main(){
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
Ith(ydot, 18) =  p[4]*y[13]*4.0+p[4]*y[46]*4.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[18]*2.0-(p[1]*y[1]*y[49]*2.0)/(p[6]+y[1])-(p[1]*y[1]*y[50]*2.0)/(p[6]+y[1])+((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 19) =  p[2]*(y[2]*y[2])*4.0+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[35]*2.0-p[2]*y[2]*y[19]*8.0+(p[1]*y[1]*y[64]*2.0)/(p[6]+y[1])+(p[1]*y[1]*y[65]*2.0)/(p[6]+y[1])+((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 20) =  p[3]*y[3]-p[3]*y[20]*2.0+p[2]*(y[2]*y[2])+p[2]*y[2]*y[51]*4.0;
Ith(ydot, 21) =  p[3]*y[3]+p[4]*y[4]-p[4]*y[21]*2.0+p[3]*y[66]*2.0;
Ith(ydot, 22) =  p[4]*y[4]+p[4]*y[5]-p[4]*y[22]*2.0+p[4]*y[80]*2.0;
Ith(ydot, 23) =  p[4]*y[5]+p[4]*y[6]-p[4]*y[23]*2.0+p[4]*y[93]*2.0;
Ith(ydot, 24) =  p[4]*y[6]+p[4]*y[7]-p[4]*y[24]*2.0+p[4]*y[105]*2.0;
Ith(ydot, 25) =  p[4]*y[7]+p[4]*y[8]-p[4]*y[25]*2.0+p[4]*y[116]*2.0;
Ith(ydot, 26) =  p[4]*y[8]+p[4]*y[9]-p[4]*y[26]*2.0+p[4]*y[126]*2.0;
Ith(ydot, 27) =  p[4]*y[9]+p[4]*y[10]-p[4]*y[27]*2.0+p[4]*y[135]*2.0;
Ith(ydot, 28) =  p[4]*y[10]+p[4]*y[11]-p[4]*y[28]*2.0+p[4]*y[143]*2.0;
Ith(ydot, 29) =  p[4]*y[11]+p[4]*y[12]-p[4]*y[29]*2.0+p[4]*y[150]*2.0;
Ith(ydot, 30) =  p[4]*y[12]+p[4]*y[13]-p[4]*y[30]*2.0+p[4]*y[156]*2.0;
Ith(ydot, 31) =  p[3]*y[3]+p[4]*y[13]+p[3]*y[76]*2.0-p[4]*y[161]*2.0;
Ith(ydot, 32) =  p[8]*y[16]+p[8]*y[168]*2.0+stm*p[5]*p[7]*y[15]-stm*p[5]*p[7]*y[32]*2.0;
Ith(ydot, 33) =  p[8]*y[16]+p[9]*y[16]-y[33]*(p[8]+p[9])*2.0+stm*p[5]*p[7]*y[15]+stm*p[5]*p[7]*y[168]*2.0;
Ith(ydot, 34) =  p[9]*y[16]+p[10]*y[17]-p[10]*y[34]*2.0+p[9]*y[170]*2.0;
Ith(ydot, 35) =  p[4]*y[61]*2.0+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[18]-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[35]-p[2]*y[2]*y[35]*4.0+(p[1]*y[1]*y[49])/(p[6]+y[1])+(p[1]*y[1]*y[50])/(p[6]+y[1])-(p[1]*y[1]*y[64])/(p[6]+y[1])-(p[1]*y[1]*y[65])/(p[6]+y[1])-((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
Ith(ydot, 36) =  -p[3]*y[36]+p[4]*y[75]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[36]+p[2]*y[2]*y[35]*2.0-(p[1]*y[1]*y[78])/(p[6]+y[1])-(p[1]*y[1]*y[79])/(p[6]+y[1]);
Ith(ydot, 37) =  p[3]*y[36]-p[4]*y[37]+p[4]*y[88]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[37]-(p[1]*y[1]*y[91])/(p[6]+y[1])-(p[1]*y[1]*y[92])/(p[6]+y[1]);
Ith(ydot, 38) =  p[4]*y[37]-p[4]*y[38]+p[4]*y[100]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[38]-(p[1]*y[1]*y[103])/(p[6]+y[1])-(p[1]*y[1]*y[104])/(p[6]+y[1]);
Ith(ydot, 39) =  p[4]*y[38]-p[4]*y[39]+p[4]*y[111]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[39]-(p[1]*y[1]*y[114])/(p[6]+y[1])-(p[1]*y[1]*y[115])/(p[6]+y[1]);
Ith(ydot, 40) =  p[4]*y[39]-p[4]*y[40]+p[4]*y[121]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[40]-(p[1]*y[1]*y[124])/(p[6]+y[1])-(p[1]*y[1]*y[125])/(p[6]+y[1]);
Ith(ydot, 41) =  p[4]*y[40]-p[4]*y[41]+p[4]*y[130]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[41]-(p[1]*y[1]*y[133])/(p[6]+y[1])-(p[1]*y[1]*y[134])/(p[6]+y[1]);
Ith(ydot, 42) =  p[4]*y[41]-p[4]*y[42]+p[4]*y[138]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[42]-(p[1]*y[1]*y[141])/(p[6]+y[1])-(p[1]*y[1]*y[142])/(p[6]+y[1]);
Ith(ydot, 43) =  p[4]*y[42]-p[4]*y[43]+p[4]*y[145]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[43]-(p[1]*y[1]*y[148])/(p[6]+y[1])-(p[1]*y[1]*y[149])/(p[6]+y[1]);
Ith(ydot, 44) =  p[4]*y[43]-p[4]*y[44]+p[4]*y[151]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[44]-(p[1]*y[1]*y[154])/(p[6]+y[1])-(p[1]*y[1]*y[155])/(p[6]+y[1]);
Ith(ydot, 45) =  p[4]*y[44]-p[4]*y[45]+p[4]*y[156]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[45]-(p[1]*y[1]*y[159])/(p[6]+y[1])-(p[1]*y[1]*y[160])/(p[6]+y[1]);
Ith(ydot, 46) =  p[4]*y[13]*-2.0+p[4]*y[30]*2.0+p[4]*y[45]-p[4]*y[46]-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[46]-(p[1]*y[1]*y[163])/(p[6]+y[1])-(p[1]*y[1]*y[164])/(p[6]+y[1]);
Ith(ydot, 47) =  p[4]*y[13]*-2.0+p[3]*y[36]-p[4]*y[46]+p[4]*y[161]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[47]-(p[1]*y[1]*y[166])/(p[6]+y[1])-(p[1]*y[1]*y[167])/(p[6]+y[1]);
Ith(ydot, 48) =  p[8]*y[49]+p[4]*y[162]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[48]-(p[1]*y[1]*y[168])/(p[6]+y[1])-(p[1]*y[1]*y[169])/(p[6]+y[1])-stm*p[5]*p[7]*y[48];
Ith(ydot, 49) =  p[4]*y[163]*2.0-y[49]*(p[8]+p[9])-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[49]-(p[1]*y[1]*y[33])/(p[6]+y[1])-(p[1]*y[1]*y[170])/(p[6]+y[1])+stm*p[5]*p[7]*y[48];
Ith(ydot, 50) =  p[9]*y[49]-p[10]*y[50]+p[4]*y[164]*2.0-(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[50]-(p[1]*y[1]*y[34])/(p[6]+y[1])-(p[1]*y[1]*y[170])/(p[6]+y[1]);
Ith(ydot, 51) =  -p[3]*y[51]-p[2]*(y[2]*y[2])*2.0+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[36]+p[2]*y[2]*y[19]*2.0-p[2]*y[2]*y[51]*4.0+(p[1]*y[1]*y[78])/(p[6]+y[1])+(p[1]*y[1]*y[79])/(p[6]+y[1]);
Ith(ydot, 52) =  p[3]*y[51]-p[4]*y[52]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[37]-p[2]*y[2]*y[52]*4.0+(p[1]*y[1]*y[91])/(p[6]+y[1])+(p[1]*y[1]*y[92])/(p[6]+y[1]);
Ith(ydot, 53) =  p[4]*y[52]-p[4]*y[53]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[38]-p[2]*y[2]*y[53]*4.0+(p[1]*y[1]*y[103])/(p[6]+y[1])+(p[1]*y[1]*y[104])/(p[6]+y[1]);
Ith(ydot, 54) =  p[4]*y[53]-p[4]*y[54]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[39]-p[2]*y[2]*y[54]*4.0+(p[1]*y[1]*y[114])/(p[6]+y[1])+(p[1]*y[1]*y[115])/(p[6]+y[1]);
Ith(ydot, 55) =  p[4]*y[54]-p[4]*y[55]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[40]-p[2]*y[2]*y[55]*4.0+(p[1]*y[1]*y[124])/(p[6]+y[1])+(p[1]*y[1]*y[125])/(p[6]+y[1]);
Ith(ydot, 56) =  p[4]*y[55]-p[4]*y[56]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[41]-p[2]*y[2]*y[56]*4.0+(p[1]*y[1]*y[133])/(p[6]+y[1])+(p[1]*y[1]*y[134])/(p[6]+y[1]);
Ith(ydot, 57) =  p[4]*y[56]-p[4]*y[57]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[42]-p[2]*y[2]*y[57]*4.0+(p[1]*y[1]*y[141])/(p[6]+y[1])+(p[1]*y[1]*y[142])/(p[6]+y[1]);
Ith(ydot, 58) =  p[4]*y[57]-p[4]*y[58]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[43]-p[2]*y[2]*y[58]*4.0+(p[1]*y[1]*y[148])/(p[6]+y[1])+(p[1]*y[1]*y[149])/(p[6]+y[1]);
Ith(ydot, 59) =  p[4]*y[58]-p[4]*y[59]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[44]-p[2]*y[2]*y[59]*4.0+(p[1]*y[1]*y[154])/(p[6]+y[1])+(p[1]*y[1]*y[155])/(p[6]+y[1]);
Ith(ydot, 60) =  p[4]*y[59]-p[4]*y[60]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[45]-p[2]*y[2]*y[60]*4.0+(p[1]*y[1]*y[159])/(p[6]+y[1])+(p[1]*y[1]*y[160])/(p[6]+y[1]);
Ith(ydot, 61) =  p[4]*y[60]-p[4]*y[61]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[46]-p[2]*y[2]*y[61]*4.0+(p[1]*y[1]*y[163])/(p[6]+y[1])+(p[1]*y[1]*y[164])/(p[6]+y[1]);
Ith(ydot, 62) =  p[3]*y[51]-p[4]*y[61]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[47]-p[2]*y[2]*y[62]*4.0+(p[1]*y[1]*y[166])/(p[6]+y[1])+(p[1]*y[1]*y[167])/(p[6]+y[1]);
Ith(ydot, 63) =  p[8]*y[64]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[48]-p[2]*y[2]*y[63]*4.0+(p[1]*y[1]*y[168])/(p[6]+y[1])+(p[1]*y[1]*y[169])/(p[6]+y[1])-stm*p[5]*p[7]*y[63];
Ith(ydot, 64) =  -y[64]*(p[8]+p[9])+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[49]-p[2]*y[2]*y[64]*4.0+(p[1]*y[1]*y[33])/(p[6]+y[1])+(p[1]*y[1]*y[170])/(p[6]+y[1])+stm*p[5]*p[7]*y[63];
Ith(ydot, 65) =  p[9]*y[64]-p[10]*y[65]+(((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1])*y[50]-p[2]*y[2]*y[65]*4.0+(p[1]*y[1]*y[34])/(p[6]+y[1])+(p[1]*y[1]*y[170])/(p[6]+y[1]);
Ith(ydot, 66) =  -p[3]*y[3]+p[3]*y[20]-p[3]*y[66]-p[4]*y[66]+p[2]*y[2]*y[52]*2.0;
Ith(ydot, 67) =  -p[3]*y[67]+p[4]*y[66]-p[4]*y[67]+p[2]*y[2]*y[53]*2.0;
Ith(ydot, 68) =  -p[3]*y[68]+p[4]*y[67]-p[4]*y[68]+p[2]*y[2]*y[54]*2.0;
Ith(ydot, 69) =  -p[3]*y[69]+p[4]*y[68]-p[4]*y[69]+p[2]*y[2]*y[55]*2.0;
Ith(ydot, 70) =  -p[3]*y[70]+p[4]*y[69]-p[4]*y[70]+p[2]*y[2]*y[56]*2.0;
Ith(ydot, 71) =  -p[3]*y[71]+p[4]*y[70]-p[4]*y[71]+p[2]*y[2]*y[57]*2.0;
Ith(ydot, 72) =  -p[3]*y[72]+p[4]*y[71]-p[4]*y[72]+p[2]*y[2]*y[58]*2.0;
Ith(ydot, 73) =  -p[3]*y[73]+p[4]*y[72]-p[4]*y[73]+p[2]*y[2]*y[59]*2.0;
Ith(ydot, 74) =  -p[3]*y[74]+p[4]*y[73]-p[4]*y[74]+p[2]*y[2]*y[60]*2.0;
Ith(ydot, 75) =  -p[3]*y[75]+p[4]*y[74]-p[4]*y[75]+p[2]*y[2]*y[61]*2.0;
Ith(ydot, 76) =  -p[3]*y[3]+p[3]*y[20]-p[3]*y[76]-p[4]*y[75]+p[2]*y[2]*y[62]*2.0;
Ith(ydot, 77) =  -p[3]*y[77]+p[8]*y[78]+p[2]*y[2]*y[63]*2.0-stm*p[5]*p[7]*y[77];
Ith(ydot, 78) =  -p[3]*y[78]-y[78]*(p[8]+p[9])+p[2]*y[2]*y[64]*2.0+stm*p[5]*p[7]*y[77];
Ith(ydot, 79) =  -p[3]*y[79]+p[9]*y[78]-p[10]*y[79]+p[2]*y[2]*y[65]*2.0;
Ith(ydot, 80) =  -p[4]*y[4]+p[4]*y[21]+p[3]*y[67]-p[4]*y[80]*2.0;
Ith(ydot, 81) =  p[3]*y[68]+p[4]*y[80]-p[4]*y[81]*2.0;
Ith(ydot, 82) =  p[3]*y[69]+p[4]*y[81]-p[4]*y[82]*2.0;
Ith(ydot, 83) =  p[3]*y[70]+p[4]*y[82]-p[4]*y[83]*2.0;
Ith(ydot, 84) =  p[3]*y[71]+p[4]*y[83]-p[4]*y[84]*2.0;
Ith(ydot, 85) =  p[3]*y[72]+p[4]*y[84]-p[4]*y[85]*2.0;
Ith(ydot, 86) =  p[3]*y[73]+p[4]*y[85]-p[4]*y[86]*2.0;
Ith(ydot, 87) =  p[3]*y[74]+p[4]*y[86]-p[4]*y[87]*2.0;
Ith(ydot, 88) =  p[3]*y[75]+p[4]*y[87]-p[4]*y[88]*2.0;
Ith(ydot, 89) =  p[3]*y[3]+p[3]*y[66]+p[3]*y[76]-p[4]*y[88]-p[4]*y[89];
Ith(ydot, 90) =  p[3]*y[77]-p[4]*y[90]+p[8]*y[91]-stm*p[5]*p[7]*y[90];
Ith(ydot, 91) =  p[3]*y[78]-p[4]*y[91]-y[91]*(p[8]+p[9])+stm*p[5]*p[7]*y[90];
Ith(ydot, 92) =  p[3]*y[79]-p[4]*y[92]+p[9]*y[91]-p[10]*y[92];
Ith(ydot, 93) =  -p[4]*y[5]+p[4]*y[22]+p[4]*y[81]-p[4]*y[93]*2.0;
Ith(ydot, 94) =  p[4]*y[82]+p[4]*y[93]-p[4]*y[94]*2.0;
Ith(ydot, 95) =  p[4]*y[83]+p[4]*y[94]-p[4]*y[95]*2.0;
Ith(ydot, 96) =  p[4]*y[84]+p[4]*y[95]-p[4]*y[96]*2.0;
Ith(ydot, 97) =  p[4]*y[85]+p[4]*y[96]-p[4]*y[97]*2.0;
Ith(ydot, 98) =  p[4]*y[86]+p[4]*y[97]-p[4]*y[98]*2.0;
Ith(ydot, 99) =  p[4]*y[87]+p[4]*y[98]-p[4]*y[99]*2.0;
Ith(ydot, 100) =  p[4]*y[88]+p[4]*y[99]-p[4]*y[100]*2.0;
Ith(ydot, 101) =  p[3]*y[67]+p[4]*y[89]-p[4]*y[100]-p[4]*y[101];
Ith(ydot, 102) =  p[4]*y[90]-p[4]*y[102]+p[8]*y[103]-stm*p[5]*p[7]*y[102];
Ith(ydot, 103) =  p[4]*y[91]-p[4]*y[103]-y[103]*(p[8]+p[9])+stm*p[5]*p[7]*y[102];
Ith(ydot, 104) =  p[4]*y[92]-p[4]*y[104]+p[9]*y[103]-p[10]*y[104];
Ith(ydot, 105) =  -p[4]*y[6]+p[4]*y[23]+p[4]*y[94]-p[4]*y[105]*2.0;
Ith(ydot, 106) =  p[4]*y[95]+p[4]*y[105]-p[4]*y[106]*2.0;
Ith(ydot, 107) =  p[4]*y[96]+p[4]*y[106]-p[4]*y[107]*2.0;
Ith(ydot, 108) =  p[4]*y[97]+p[4]*y[107]-p[4]*y[108]*2.0;
Ith(ydot, 109) =  p[4]*y[98]+p[4]*y[108]-p[4]*y[109]*2.0;
Ith(ydot, 110) =  p[4]*y[99]+p[4]*y[109]-p[4]*y[110]*2.0;
Ith(ydot, 111) =  p[4]*y[100]+p[4]*y[110]-p[4]*y[111]*2.0;
Ith(ydot, 112) =  p[3]*y[68]+p[4]*y[101]-p[4]*y[111]-p[4]*y[112];
Ith(ydot, 113) =  p[4]*y[102]-p[4]*y[113]+p[8]*y[114]-stm*p[5]*p[7]*y[113];
Ith(ydot, 114) =  p[4]*y[103]-p[4]*y[114]-y[114]*(p[8]+p[9])+stm*p[5]*p[7]*y[113];
Ith(ydot, 115) =  p[4]*y[104]-p[4]*y[115]+p[9]*y[114]-p[10]*y[115];
Ith(ydot, 116) =  -p[4]*y[7]+p[4]*y[24]+p[4]*y[106]-p[4]*y[116]*2.0;
Ith(ydot, 117) =  p[4]*y[107]+p[4]*y[116]-p[4]*y[117]*2.0;
Ith(ydot, 118) =  p[4]*y[108]+p[4]*y[117]-p[4]*y[118]*2.0;
Ith(ydot, 119) =  p[4]*y[109]+p[4]*y[118]-p[4]*y[119]*2.0;
Ith(ydot, 120) =  p[4]*y[110]+p[4]*y[119]-p[4]*y[120]*2.0;
Ith(ydot, 121) =  p[4]*y[111]+p[4]*y[120]-p[4]*y[121]*2.0;
Ith(ydot, 122) =  p[3]*y[69]+p[4]*y[112]-p[4]*y[121]-p[4]*y[122];
Ith(ydot, 123) =  p[4]*y[113]-p[4]*y[123]+p[8]*y[124]-stm*p[5]*p[7]*y[123];
Ith(ydot, 124) =  p[4]*y[114]-p[4]*y[124]-y[124]*(p[8]+p[9])+stm*p[5]*p[7]*y[123];
Ith(ydot, 125) =  p[4]*y[115]-p[4]*y[125]+p[9]*y[124]-p[10]*y[125];
Ith(ydot, 126) =  -p[4]*y[8]+p[4]*y[25]+p[4]*y[117]-p[4]*y[126]*2.0;
Ith(ydot, 127) =  p[4]*y[118]+p[4]*y[126]-p[4]*y[127]*2.0;
Ith(ydot, 128) =  p[4]*y[119]+p[4]*y[127]-p[4]*y[128]*2.0;
Ith(ydot, 129) =  p[4]*y[120]+p[4]*y[128]-p[4]*y[129]*2.0;
Ith(ydot, 130) =  p[4]*y[121]+p[4]*y[129]-p[4]*y[130]*2.0;
Ith(ydot, 131) =  p[3]*y[70]+p[4]*y[122]-p[4]*y[130]-p[4]*y[131];
Ith(ydot, 132) =  p[4]*y[123]-p[4]*y[132]+p[8]*y[133]-stm*p[5]*p[7]*y[132];
Ith(ydot, 133) =  p[4]*y[124]-p[4]*y[133]-y[133]*(p[8]+p[9])+stm*p[5]*p[7]*y[132];
Ith(ydot, 134) =  p[4]*y[125]-p[4]*y[134]+p[9]*y[133]-p[10]*y[134];
Ith(ydot, 135) =  -p[4]*y[9]+p[4]*y[26]+p[4]*y[127]-p[4]*y[135]*2.0;
Ith(ydot, 136) =  p[4]*y[128]+p[4]*y[135]-p[4]*y[136]*2.0;
Ith(ydot, 137) =  p[4]*y[129]+p[4]*y[136]-p[4]*y[137]*2.0;
Ith(ydot, 138) =  p[4]*y[130]+p[4]*y[137]-p[4]*y[138]*2.0;
Ith(ydot, 139) =  p[3]*y[71]+p[4]*y[131]-p[4]*y[138]-p[4]*y[139];
Ith(ydot, 140) =  p[4]*y[132]-p[4]*y[140]+p[8]*y[141]-stm*p[5]*p[7]*y[140];
Ith(ydot, 141) =  p[4]*y[133]-p[4]*y[141]-y[141]*(p[8]+p[9])+stm*p[5]*p[7]*y[140];
Ith(ydot, 142) =  p[4]*y[134]-p[4]*y[142]+p[9]*y[141]-p[10]*y[142];
Ith(ydot, 143) =  -p[4]*y[10]+p[4]*y[27]+p[4]*y[136]-p[4]*y[143]*2.0;
Ith(ydot, 144) =  p[4]*y[137]+p[4]*y[143]-p[4]*y[144]*2.0;
Ith(ydot, 145) =  p[4]*y[138]+p[4]*y[144]-p[4]*y[145]*2.0;
Ith(ydot, 146) =  p[3]*y[72]+p[4]*y[139]-p[4]*y[145]-p[4]*y[146];
Ith(ydot, 147) =  p[4]*y[140]-p[4]*y[147]+p[8]*y[148]-stm*p[5]*p[7]*y[147];
Ith(ydot, 148) =  p[4]*y[141]-p[4]*y[148]-y[148]*(p[8]+p[9])+stm*p[5]*p[7]*y[147];
Ith(ydot, 149) =  p[4]*y[142]-p[4]*y[149]+p[9]*y[148]-p[10]*y[149];
Ith(ydot, 150) =  -p[4]*y[11]+p[4]*y[28]+p[4]*y[144]-p[4]*y[150]*2.0;
Ith(ydot, 151) =  p[4]*y[145]+p[4]*y[150]-p[4]*y[151]*2.0;
Ith(ydot, 152) =  p[3]*y[73]+p[4]*y[146]-p[4]*y[151]-p[4]*y[152];
Ith(ydot, 153) =  p[4]*y[147]-p[4]*y[153]+p[8]*y[154]-stm*p[5]*p[7]*y[153];
Ith(ydot, 154) =  p[4]*y[148]-p[4]*y[154]-y[154]*(p[8]+p[9])+stm*p[5]*p[7]*y[153];
Ith(ydot, 155) =  p[4]*y[149]-p[4]*y[155]+p[9]*y[154]-p[10]*y[155];
Ith(ydot, 156) =  -p[4]*y[12]+p[4]*y[29]+p[4]*y[151]-p[4]*y[156]*2.0;
Ith(ydot, 157) =  p[3]*y[74]+p[4]*y[152]-p[4]*y[156]-p[4]*y[157];
Ith(ydot, 158) =  p[4]*y[153]-p[4]*y[158]+p[8]*y[159]-stm*p[5]*p[7]*y[158];
Ith(ydot, 159) =  p[4]*y[154]-p[4]*y[159]-y[159]*(p[8]+p[9])+stm*p[5]*p[7]*y[158];
Ith(ydot, 160) =  p[4]*y[155]-p[4]*y[160]+p[9]*y[159]-p[10]*y[160];
Ith(ydot, 161) =  p[4]*y[13]-p[4]*y[30]+p[3]*y[75]+p[4]*y[157]-p[4]*y[161];
Ith(ydot, 162) =  p[4]*y[158]-p[4]*y[162]+p[8]*y[163]-stm*p[5]*p[7]*y[162];
Ith(ydot, 163) =  p[4]*y[159]-p[4]*y[163]-y[163]*(p[8]+p[9])+stm*p[5]*p[7]*y[162];
Ith(ydot, 164) =  p[4]*y[160]-p[4]*y[164]+p[9]*y[163]-p[10]*y[164];
Ith(ydot, 165) =  p[3]*y[77]-p[4]*y[162]+p[8]*y[166]-stm*p[5]*p[7]*y[165];
Ith(ydot, 166) =  p[3]*y[78]-p[4]*y[163]-y[166]*(p[8]+p[9])+stm*p[5]*p[7]*y[165];
Ith(ydot, 167) =  p[3]*y[79]-p[4]*y[164]+p[9]*y[166]-p[10]*y[167];
Ith(ydot, 168) =  -p[8]*y[16]+p[8]*y[33]-y[168]*(p[8]+p[9])-stm*p[5]*p[7]*y[15]+stm*p[5]*p[7]*y[32]-stm*p[5]*p[7]*y[168];
Ith(ydot, 169) =  p[9]*y[168]+p[8]*y[170]-p[10]*y[169]-stm*p[5]*p[7]*y[169];
Ith(ydot, 170) =  -p[9]*y[16]+p[9]*y[33]-p[10]*y[170]-y[170]*(p[8]+p[9])+stm*p[5]*p[7]*y[169];


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
IJth(J, 18, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[18]*2.0-(p[1]*y[49]*2.0)/(p[6]+y[1])-(p[1]*y[50]*2.0)/(p[6]+y[1])+((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[49]*2.0+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[50]*2.0-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 18, 13) =  p[4]*4.0;
IJth(J, 18, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[18]*-2.0+(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 18, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[18]*-2.0+(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 18, 18) =  ((y[16]+y[17])*p[1]*-2.0)/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1]*2.0;
IJth(J, 18, 46) =  p[4]*4.0;
IJth(J, 18, 49) =  (p[1]*y[1]*-2.0)/(p[6]+y[1]);
IJth(J, 18, 50) =  (p[1]*y[1]*-2.0)/(p[6]+y[1]);
IJth(J, 19, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[35]*-2.0+(p[1]*y[64]*2.0)/(p[6]+y[1])+(p[1]*y[65]*2.0)/(p[6]+y[1])+((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[64]*2.0-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[65]*2.0-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 19, 2) =  p[2]*y[2]*8.0-p[2]*y[19]*8.0;
IJth(J, 19, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[35]*2.0+(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 19, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[35]*2.0+(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 19, 19) =  p[2]*y[2]*-8.0;
IJth(J, 19, 35) =  ((y[16]+y[17])*p[1]*2.0)/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1]*2.0;
IJth(J, 19, 64) =  (p[1]*y[1]*2.0)/(p[6]+y[1]);
IJth(J, 19, 65) =  (p[1]*y[1]*2.0)/(p[6]+y[1]);
IJth(J, 20, 2) =  p[2]*y[2]*2.0+p[2]*y[51]*4.0;
IJth(J, 20, 3) =  p[3];
IJth(J, 20, 20) =  p[3]*-2.0;
IJth(J, 20, 51) =  p[2]*y[2]*4.0;
IJth(J, 21, 3) =  p[3];
IJth(J, 21, 4) =  p[4];
IJth(J, 21, 21) =  p[4]*-2.0;
IJth(J, 21, 66) =  p[3]*2.0;
IJth(J, 22, 4) =  p[4];
IJth(J, 22, 5) =  p[4];
IJth(J, 22, 22) =  p[4]*-2.0;
IJth(J, 22, 80) =  p[4]*2.0;
IJth(J, 23, 5) =  p[4];
IJth(J, 23, 6) =  p[4];
IJth(J, 23, 23) =  p[4]*-2.0;
IJth(J, 23, 93) =  p[4]*2.0;
IJth(J, 24, 6) =  p[4];
IJth(J, 24, 7) =  p[4];
IJth(J, 24, 24) =  p[4]*-2.0;
IJth(J, 24, 105) =  p[4]*2.0;
IJth(J, 25, 7) =  p[4];
IJth(J, 25, 8) =  p[4];
IJth(J, 25, 25) =  p[4]*-2.0;
IJth(J, 25, 116) =  p[4]*2.0;
IJth(J, 26, 8) =  p[4];
IJth(J, 26, 9) =  p[4];
IJth(J, 26, 26) =  p[4]*-2.0;
IJth(J, 26, 126) =  p[4]*2.0;
IJth(J, 27, 9) =  p[4];
IJth(J, 27, 10) =  p[4];
IJth(J, 27, 27) =  p[4]*-2.0;
IJth(J, 27, 135) =  p[4]*2.0;
IJth(J, 28, 10) =  p[4];
IJth(J, 28, 11) =  p[4];
IJth(J, 28, 28) =  p[4]*-2.0;
IJth(J, 28, 143) =  p[4]*2.0;
IJth(J, 29, 11) =  p[4];
IJth(J, 29, 12) =  p[4];
IJth(J, 29, 29) =  p[4]*-2.0;
IJth(J, 29, 150) =  p[4]*2.0;
IJth(J, 30, 12) =  p[4];
IJth(J, 30, 13) =  p[4];
IJth(J, 30, 30) =  p[4]*-2.0;
IJth(J, 30, 156) =  p[4]*2.0;
IJth(J, 31, 3) =  p[3];
IJth(J, 31, 13) =  p[4];
IJth(J, 31, 76) =  p[3]*2.0;
IJth(J, 31, 161) =  p[4]*-2.0;
IJth(J, 32, 15) =  stm*p[5]*p[7];
IJth(J, 32, 16) =  p[8];
IJth(J, 32, 32) =  stm*p[5]*p[7]*-2.0;
IJth(J, 32, 168) =  p[8]*2.0;
IJth(J, 33, 15) =  stm*p[5]*p[7];
IJth(J, 33, 16) =  p[8]+p[9];
IJth(J, 33, 33) =  p[8]*-2.0-p[9]*2.0;
IJth(J, 33, 168) =  stm*p[5]*p[7]*2.0;
IJth(J, 34, 16) =  p[9];
IJth(J, 34, 17) =  p[10];
IJth(J, 34, 34) =  p[10]*-2.0;
IJth(J, 34, 170) =  p[9]*2.0;
IJth(J, 35, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[18]+(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[35]+(p[1]*y[49])/(p[6]+y[1])+(p[1]*y[50])/(p[6]+y[1])-(p[1]*y[64])/(p[6]+y[1])-(p[1]*y[65])/(p[6]+y[1])-((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[49]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[50]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[64]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[65]+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 35, 2) =  p[2]*y[35]*-4.0;
IJth(J, 35, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[18]-(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[35]-(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 35, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[18]-(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[35]-(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 35, 18) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 35, 35) =  p[2]*y[2]*-4.0-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 35, 49) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 35, 50) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 35, 61) =  p[4]*2.0;
IJth(J, 35, 64) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 35, 65) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 36, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[36]-(p[1]*y[78])/(p[6]+y[1])-(p[1]*y[79])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[78]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[79];
IJth(J, 36, 2) =  p[2]*y[35]*2.0;
IJth(J, 36, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[36];
IJth(J, 36, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[36];
IJth(J, 36, 35) =  p[2]*y[2]*2.0;
IJth(J, 36, 36) =  -p[3]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 36, 75) =  p[4]*2.0;
IJth(J, 36, 78) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 36, 79) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 37, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[37]-(p[1]*y[91])/(p[6]+y[1])-(p[1]*y[92])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[91]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[92];
IJth(J, 37, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[37];
IJth(J, 37, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[37];
IJth(J, 37, 36) =  p[3];
IJth(J, 37, 37) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 37, 88) =  p[4]*2.0;
IJth(J, 37, 91) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 37, 92) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 38, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[38]-(p[1]*y[103])/(p[6]+y[1])-(p[1]*y[104])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[103]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[104];
IJth(J, 38, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[38];
IJth(J, 38, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[38];
IJth(J, 38, 37) =  p[4];
IJth(J, 38, 38) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 38, 100) =  p[4]*2.0;
IJth(J, 38, 103) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 38, 104) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 39, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[39]-(p[1]*y[114])/(p[6]+y[1])-(p[1]*y[115])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[114]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[115];
IJth(J, 39, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[39];
IJth(J, 39, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[39];
IJth(J, 39, 38) =  p[4];
IJth(J, 39, 39) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 39, 111) =  p[4]*2.0;
IJth(J, 39, 114) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 39, 115) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 40, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[40]-(p[1]*y[124])/(p[6]+y[1])-(p[1]*y[125])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[124]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[125];
IJth(J, 40, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[40];
IJth(J, 40, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[40];
IJth(J, 40, 39) =  p[4];
IJth(J, 40, 40) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 40, 121) =  p[4]*2.0;
IJth(J, 40, 124) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 40, 125) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 41, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[41]-(p[1]*y[133])/(p[6]+y[1])-(p[1]*y[134])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[133]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[134];
IJth(J, 41, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[41];
IJth(J, 41, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[41];
IJth(J, 41, 40) =  p[4];
IJth(J, 41, 41) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 41, 130) =  p[4]*2.0;
IJth(J, 41, 133) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 41, 134) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 42, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[42]-(p[1]*y[141])/(p[6]+y[1])-(p[1]*y[142])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[141]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[142];
IJth(J, 42, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[42];
IJth(J, 42, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[42];
IJth(J, 42, 41) =  p[4];
IJth(J, 42, 42) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 42, 138) =  p[4]*2.0;
IJth(J, 42, 141) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 42, 142) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 43, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[43]-(p[1]*y[148])/(p[6]+y[1])-(p[1]*y[149])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[148]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[149];
IJth(J, 43, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[43];
IJth(J, 43, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[43];
IJth(J, 43, 42) =  p[4];
IJth(J, 43, 43) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 43, 145) =  p[4]*2.0;
IJth(J, 43, 148) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 43, 149) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 44, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[44]-(p[1]*y[154])/(p[6]+y[1])-(p[1]*y[155])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[154]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[155];
IJth(J, 44, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[44];
IJth(J, 44, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[44];
IJth(J, 44, 43) =  p[4];
IJth(J, 44, 44) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 44, 151) =  p[4]*2.0;
IJth(J, 44, 154) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 44, 155) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 45, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[45]-(p[1]*y[159])/(p[6]+y[1])-(p[1]*y[160])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[159]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[160];
IJth(J, 45, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[45];
IJth(J, 45, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[45];
IJth(J, 45, 44) =  p[4];
IJth(J, 45, 45) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 45, 156) =  p[4]*2.0;
IJth(J, 45, 159) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 45, 160) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 46, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[46]-(p[1]*y[163])/(p[6]+y[1])-(p[1]*y[164])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[163]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[164];
IJth(J, 46, 13) =  p[4]*-2.0;
IJth(J, 46, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[46];
IJth(J, 46, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[46];
IJth(J, 46, 30) =  p[4]*2.0;
IJth(J, 46, 45) =  p[4];
IJth(J, 46, 46) =  -p[4]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 46, 163) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 46, 164) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 47, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[47]-(p[1]*y[166])/(p[6]+y[1])-(p[1]*y[167])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[166]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[167];
IJth(J, 47, 13) =  p[4]*-2.0;
IJth(J, 47, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[47];
IJth(J, 47, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[47];
IJth(J, 47, 36) =  p[3];
IJth(J, 47, 46) =  -p[4];
IJth(J, 47, 47) =  -((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 47, 161) =  p[4]*2.0;
IJth(J, 47, 166) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 47, 167) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 48, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[48]-(p[1]*y[168])/(p[6]+y[1])-(p[1]*y[169])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[168]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[169];
IJth(J, 48, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[48];
IJth(J, 48, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[48];
IJth(J, 48, 48) =  -stm*p[5]*p[7]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 48, 49) =  p[8];
IJth(J, 48, 162) =  p[4]*2.0;
IJth(J, 48, 168) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 48, 169) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 49, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[49]-(p[1]*y[33])/(p[6]+y[1])-(p[1]*y[170])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[33]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[170];
IJth(J, 49, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[49];
IJth(J, 49, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[49];
IJth(J, 49, 33) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 49, 48) =  stm*p[5]*p[7];
IJth(J, 49, 49) =  -p[8]-p[9]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 49, 163) =  p[4]*2.0;
IJth(J, 49, 170) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 50, 1) =  (1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[50]-(p[1]*y[34])/(p[6]+y[1])-(p[1]*y[170])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[34]+1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[170];
IJth(J, 50, 16) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[50];
IJth(J, 50, 17) =  -(p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[50];
IJth(J, 50, 34) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 50, 49) =  p[9];
IJth(J, 50, 50) =  -p[10]-((y[16]+y[17])*p[1])/(p[6]+y[1])+1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 50, 164) =  p[4]*2.0;
IJth(J, 50, 170) =  -(p[1]*y[1])/(p[6]+y[1]);
IJth(J, 51, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[36]+(p[1]*y[78])/(p[6]+y[1])+(p[1]*y[79])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[78]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[79];
IJth(J, 51, 2) =  p[2]*y[2]*-4.0+p[2]*y[19]*2.0-p[2]*y[51]*4.0;
IJth(J, 51, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[36];
IJth(J, 51, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[36];
IJth(J, 51, 19) =  p[2]*y[2]*2.0;
IJth(J, 51, 36) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 51, 51) =  p[2]*y[2]*-4.0-p[3];
IJth(J, 51, 78) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 51, 79) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 52, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[37]+(p[1]*y[91])/(p[6]+y[1])+(p[1]*y[92])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[91]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[92];
IJth(J, 52, 2) =  p[2]*y[52]*-4.0;
IJth(J, 52, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[37];
IJth(J, 52, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[37];
IJth(J, 52, 37) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 52, 51) =  p[3];
IJth(J, 52, 52) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 52, 91) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 52, 92) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 53, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[38]+(p[1]*y[103])/(p[6]+y[1])+(p[1]*y[104])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[103]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[104];
IJth(J, 53, 2) =  p[2]*y[53]*-4.0;
IJth(J, 53, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[38];
IJth(J, 53, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[38];
IJth(J, 53, 38) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 53, 52) =  p[4];
IJth(J, 53, 53) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 53, 103) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 53, 104) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 54, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[39]+(p[1]*y[114])/(p[6]+y[1])+(p[1]*y[115])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[114]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[115];
IJth(J, 54, 2) =  p[2]*y[54]*-4.0;
IJth(J, 54, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[39];
IJth(J, 54, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[39];
IJth(J, 54, 39) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 54, 53) =  p[4];
IJth(J, 54, 54) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 54, 114) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 54, 115) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 55, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[40]+(p[1]*y[124])/(p[6]+y[1])+(p[1]*y[125])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[124]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[125];
IJth(J, 55, 2) =  p[2]*y[55]*-4.0;
IJth(J, 55, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[40];
IJth(J, 55, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[40];
IJth(J, 55, 40) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 55, 54) =  p[4];
IJth(J, 55, 55) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 55, 124) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 55, 125) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 56, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[41]+(p[1]*y[133])/(p[6]+y[1])+(p[1]*y[134])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[133]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[134];
IJth(J, 56, 2) =  p[2]*y[56]*-4.0;
IJth(J, 56, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[41];
IJth(J, 56, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[41];
IJth(J, 56, 41) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 56, 55) =  p[4];
IJth(J, 56, 56) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 56, 133) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 56, 134) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 57, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[42]+(p[1]*y[141])/(p[6]+y[1])+(p[1]*y[142])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[141]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[142];
IJth(J, 57, 2) =  p[2]*y[57]*-4.0;
IJth(J, 57, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[42];
IJth(J, 57, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[42];
IJth(J, 57, 42) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 57, 56) =  p[4];
IJth(J, 57, 57) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 57, 141) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 57, 142) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 58, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[43]+(p[1]*y[148])/(p[6]+y[1])+(p[1]*y[149])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[148]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[149];
IJth(J, 58, 2) =  p[2]*y[58]*-4.0;
IJth(J, 58, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[43];
IJth(J, 58, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[43];
IJth(J, 58, 43) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 58, 57) =  p[4];
IJth(J, 58, 58) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 58, 148) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 58, 149) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 59, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[44]+(p[1]*y[154])/(p[6]+y[1])+(p[1]*y[155])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[154]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[155];
IJth(J, 59, 2) =  p[2]*y[59]*-4.0;
IJth(J, 59, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[44];
IJth(J, 59, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[44];
IJth(J, 59, 44) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 59, 58) =  p[4];
IJth(J, 59, 59) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 59, 154) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 59, 155) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 60, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[45]+(p[1]*y[159])/(p[6]+y[1])+(p[1]*y[160])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[159]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[160];
IJth(J, 60, 2) =  p[2]*y[60]*-4.0;
IJth(J, 60, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[45];
IJth(J, 60, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[45];
IJth(J, 60, 45) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 60, 59) =  p[4];
IJth(J, 60, 60) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 60, 159) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 60, 160) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 61, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[46]+(p[1]*y[163])/(p[6]+y[1])+(p[1]*y[164])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[163]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[164];
IJth(J, 61, 2) =  p[2]*y[61]*-4.0;
IJth(J, 61, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[46];
IJth(J, 61, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[46];
IJth(J, 61, 46) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 61, 60) =  p[4];
IJth(J, 61, 61) =  p[2]*y[2]*-4.0-p[4];
IJth(J, 61, 163) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 61, 164) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 62, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[47]+(p[1]*y[166])/(p[6]+y[1])+(p[1]*y[167])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[166]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[167];
IJth(J, 62, 2) =  p[2]*y[62]*-4.0;
IJth(J, 62, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[47];
IJth(J, 62, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[47];
IJth(J, 62, 47) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 62, 51) =  p[3];
IJth(J, 62, 61) =  -p[4];
IJth(J, 62, 62) =  p[2]*y[2]*-4.0;
IJth(J, 62, 166) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 62, 167) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 63, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[48]+(p[1]*y[168])/(p[6]+y[1])+(p[1]*y[169])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[168]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[169];
IJth(J, 63, 2) =  p[2]*y[63]*-4.0;
IJth(J, 63, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[48];
IJth(J, 63, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[48];
IJth(J, 63, 48) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 63, 63) =  p[2]*y[2]*-4.0-stm*p[5]*p[7];
IJth(J, 63, 64) =  p[8];
IJth(J, 63, 168) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 63, 169) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 64, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[49]+(p[1]*y[33])/(p[6]+y[1])+(p[1]*y[170])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[33]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[170];
IJth(J, 64, 2) =  p[2]*y[64]*-4.0;
IJth(J, 64, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[49];
IJth(J, 64, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[49];
IJth(J, 64, 33) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 64, 49) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 64, 63) =  stm*p[5]*p[7];
IJth(J, 64, 64) =  p[2]*y[2]*-4.0-p[8]-p[9];
IJth(J, 64, 170) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 65, 1) =  -(1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*2.0-1.0/pow(p[6]+y[1],3.0)*(y[16]+y[17])*p[1]*y[1]*2.0)*y[50]+(p[1]*y[34])/(p[6]+y[1])+(p[1]*y[170])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[34]-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1]*y[170];
IJth(J, 65, 2) =  p[2]*y[65]*-4.0;
IJth(J, 65, 16) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[50];
IJth(J, 65, 17) =  (p[1]/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*p[1]*y[1])*y[50];
IJth(J, 65, 34) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 65, 50) =  ((y[16]+y[17])*p[1])/(p[6]+y[1])-1.0/pow(p[6]+y[1],2.0)*(y[16]+y[17])*p[1]*y[1];
IJth(J, 65, 64) =  p[9];
IJth(J, 65, 65) =  p[2]*y[2]*-4.0-p[10];
IJth(J, 65, 170) =  (p[1]*y[1])/(p[6]+y[1]);
IJth(J, 66, 2) =  p[2]*y[52]*2.0;
IJth(J, 66, 3) =  -p[3];
IJth(J, 66, 20) =  p[3];
IJth(J, 66, 52) =  p[2]*y[2]*2.0;
IJth(J, 66, 66) =  -p[3]-p[4];
IJth(J, 67, 2) =  p[2]*y[53]*2.0;
IJth(J, 67, 53) =  p[2]*y[2]*2.0;
IJth(J, 67, 66) =  p[4];
IJth(J, 67, 67) =  -p[3]-p[4];
IJth(J, 68, 2) =  p[2]*y[54]*2.0;
IJth(J, 68, 54) =  p[2]*y[2]*2.0;
IJth(J, 68, 67) =  p[4];
IJth(J, 68, 68) =  -p[3]-p[4];
IJth(J, 69, 2) =  p[2]*y[55]*2.0;
IJth(J, 69, 55) =  p[2]*y[2]*2.0;
IJth(J, 69, 68) =  p[4];
IJth(J, 69, 69) =  -p[3]-p[4];
IJth(J, 70, 2) =  p[2]*y[56]*2.0;
IJth(J, 70, 56) =  p[2]*y[2]*2.0;
IJth(J, 70, 69) =  p[4];
IJth(J, 70, 70) =  -p[3]-p[4];
IJth(J, 71, 2) =  p[2]*y[57]*2.0;
IJth(J, 71, 57) =  p[2]*y[2]*2.0;
IJth(J, 71, 70) =  p[4];
IJth(J, 71, 71) =  -p[3]-p[4];
IJth(J, 72, 2) =  p[2]*y[58]*2.0;
IJth(J, 72, 58) =  p[2]*y[2]*2.0;
IJth(J, 72, 71) =  p[4];
IJth(J, 72, 72) =  -p[3]-p[4];
IJth(J, 73, 2) =  p[2]*y[59]*2.0;
IJth(J, 73, 59) =  p[2]*y[2]*2.0;
IJth(J, 73, 72) =  p[4];
IJth(J, 73, 73) =  -p[3]-p[4];
IJth(J, 74, 2) =  p[2]*y[60]*2.0;
IJth(J, 74, 60) =  p[2]*y[2]*2.0;
IJth(J, 74, 73) =  p[4];
IJth(J, 74, 74) =  -p[3]-p[4];
IJth(J, 75, 2) =  p[2]*y[61]*2.0;
IJth(J, 75, 61) =  p[2]*y[2]*2.0;
IJth(J, 75, 74) =  p[4];
IJth(J, 75, 75) =  -p[3]-p[4];
IJth(J, 76, 2) =  p[2]*y[62]*2.0;
IJth(J, 76, 3) =  -p[3];
IJth(J, 76, 20) =  p[3];
IJth(J, 76, 62) =  p[2]*y[2]*2.0;
IJth(J, 76, 75) =  -p[4];
IJth(J, 76, 76) =  -p[3];
IJth(J, 77, 2) =  p[2]*y[63]*2.0;
IJth(J, 77, 63) =  p[2]*y[2]*2.0;
IJth(J, 77, 77) =  -p[3]-stm*p[5]*p[7];
IJth(J, 77, 78) =  p[8];
IJth(J, 78, 2) =  p[2]*y[64]*2.0;
IJth(J, 78, 64) =  p[2]*y[2]*2.0;
IJth(J, 78, 77) =  stm*p[5]*p[7];
IJth(J, 78, 78) =  -p[3]-p[8]-p[9];
IJth(J, 79, 2) =  p[2]*y[65]*2.0;
IJth(J, 79, 65) =  p[2]*y[2]*2.0;
IJth(J, 79, 78) =  p[9];
IJth(J, 79, 79) =  -p[3]-p[10];
IJth(J, 80, 4) =  -p[4];
IJth(J, 80, 21) =  p[4];
IJth(J, 80, 67) =  p[3];
IJth(J, 80, 80) =  p[4]*-2.0;
IJth(J, 81, 68) =  p[3];
IJth(J, 81, 80) =  p[4];
IJth(J, 81, 81) =  p[4]*-2.0;
IJth(J, 82, 69) =  p[3];
IJth(J, 82, 81) =  p[4];
IJth(J, 82, 82) =  p[4]*-2.0;
IJth(J, 83, 70) =  p[3];
IJth(J, 83, 82) =  p[4];
IJth(J, 83, 83) =  p[4]*-2.0;
IJth(J, 84, 71) =  p[3];
IJth(J, 84, 83) =  p[4];
IJth(J, 84, 84) =  p[4]*-2.0;
IJth(J, 85, 72) =  p[3];
IJth(J, 85, 84) =  p[4];
IJth(J, 85, 85) =  p[4]*-2.0;
IJth(J, 86, 73) =  p[3];
IJth(J, 86, 85) =  p[4];
IJth(J, 86, 86) =  p[4]*-2.0;
IJth(J, 87, 74) =  p[3];
IJth(J, 87, 86) =  p[4];
IJth(J, 87, 87) =  p[4]*-2.0;
IJth(J, 88, 75) =  p[3];
IJth(J, 88, 87) =  p[4];
IJth(J, 88, 88) =  p[4]*-2.0;
IJth(J, 89, 3) =  p[3];
IJth(J, 89, 66) =  p[3];
IJth(J, 89, 76) =  p[3];
IJth(J, 89, 88) =  -p[4];
IJth(J, 89, 89) =  -p[4];
IJth(J, 90, 77) =  p[3];
IJth(J, 90, 90) =  -p[4]-stm*p[5]*p[7];
IJth(J, 90, 91) =  p[8];
IJth(J, 91, 78) =  p[3];
IJth(J, 91, 90) =  stm*p[5]*p[7];
IJth(J, 91, 91) =  -p[4]-p[8]-p[9];
IJth(J, 92, 79) =  p[3];
IJth(J, 92, 91) =  p[9];
IJth(J, 92, 92) =  -p[4]-p[10];
IJth(J, 93, 5) =  -p[4];
IJth(J, 93, 22) =  p[4];
IJth(J, 93, 81) =  p[4];
IJth(J, 93, 93) =  p[4]*-2.0;
IJth(J, 94, 82) =  p[4];
IJth(J, 94, 93) =  p[4];
IJth(J, 94, 94) =  p[4]*-2.0;
IJth(J, 95, 83) =  p[4];
IJth(J, 95, 94) =  p[4];
IJth(J, 95, 95) =  p[4]*-2.0;
IJth(J, 96, 84) =  p[4];
IJth(J, 96, 95) =  p[4];
IJth(J, 96, 96) =  p[4]*-2.0;
IJth(J, 97, 85) =  p[4];
IJth(J, 97, 96) =  p[4];
IJth(J, 97, 97) =  p[4]*-2.0;
IJth(J, 98, 86) =  p[4];
IJth(J, 98, 97) =  p[4];
IJth(J, 98, 98) =  p[4]*-2.0;
IJth(J, 99, 87) =  p[4];
IJth(J, 99, 98) =  p[4];
IJth(J, 99, 99) =  p[4]*-2.0;
IJth(J, 100, 88) =  p[4];
IJth(J, 100, 99) =  p[4];
IJth(J, 100, 100) =  p[4]*-2.0;
IJth(J, 101, 67) =  p[3];
IJth(J, 101, 89) =  p[4];
IJth(J, 101, 100) =  -p[4];
IJth(J, 101, 101) =  -p[4];
IJth(J, 102, 90) =  p[4];
IJth(J, 102, 102) =  -p[4]-stm*p[5]*p[7];
IJth(J, 102, 103) =  p[8];
IJth(J, 103, 91) =  p[4];
IJth(J, 103, 102) =  stm*p[5]*p[7];
IJth(J, 103, 103) =  -p[4]-p[8]-p[9];
IJth(J, 104, 92) =  p[4];
IJth(J, 104, 103) =  p[9];
IJth(J, 104, 104) =  -p[4]-p[10];
IJth(J, 105, 6) =  -p[4];
IJth(J, 105, 23) =  p[4];
IJth(J, 105, 94) =  p[4];
IJth(J, 105, 105) =  p[4]*-2.0;
IJth(J, 106, 95) =  p[4];
IJth(J, 106, 105) =  p[4];
IJth(J, 106, 106) =  p[4]*-2.0;
IJth(J, 107, 96) =  p[4];
IJth(J, 107, 106) =  p[4];
IJth(J, 107, 107) =  p[4]*-2.0;
IJth(J, 108, 97) =  p[4];
IJth(J, 108, 107) =  p[4];
IJth(J, 108, 108) =  p[4]*-2.0;
IJth(J, 109, 98) =  p[4];
IJth(J, 109, 108) =  p[4];
IJth(J, 109, 109) =  p[4]*-2.0;
IJth(J, 110, 99) =  p[4];
IJth(J, 110, 109) =  p[4];
IJth(J, 110, 110) =  p[4]*-2.0;
IJth(J, 111, 100) =  p[4];
IJth(J, 111, 110) =  p[4];
IJth(J, 111, 111) =  p[4]*-2.0;
IJth(J, 112, 68) =  p[3];
IJth(J, 112, 101) =  p[4];
IJth(J, 112, 111) =  -p[4];
IJth(J, 112, 112) =  -p[4];
IJth(J, 113, 102) =  p[4];
IJth(J, 113, 113) =  -p[4]-stm*p[5]*p[7];
IJth(J, 113, 114) =  p[8];
IJth(J, 114, 103) =  p[4];
IJth(J, 114, 113) =  stm*p[5]*p[7];
IJth(J, 114, 114) =  -p[4]-p[8]-p[9];
IJth(J, 115, 104) =  p[4];
IJth(J, 115, 114) =  p[9];
IJth(J, 115, 115) =  -p[4]-p[10];
IJth(J, 116, 7) =  -p[4];
IJth(J, 116, 24) =  p[4];
IJth(J, 116, 106) =  p[4];
IJth(J, 116, 116) =  p[4]*-2.0;
IJth(J, 117, 107) =  p[4];
IJth(J, 117, 116) =  p[4];
IJth(J, 117, 117) =  p[4]*-2.0;
IJth(J, 118, 108) =  p[4];
IJth(J, 118, 117) =  p[4];
IJth(J, 118, 118) =  p[4]*-2.0;
IJth(J, 119, 109) =  p[4];
IJth(J, 119, 118) =  p[4];
IJth(J, 119, 119) =  p[4]*-2.0;
IJth(J, 120, 110) =  p[4];
IJth(J, 120, 119) =  p[4];
IJth(J, 120, 120) =  p[4]*-2.0;
IJth(J, 121, 111) =  p[4];
IJth(J, 121, 120) =  p[4];
IJth(J, 121, 121) =  p[4]*-2.0;
IJth(J, 122, 69) =  p[3];
IJth(J, 122, 112) =  p[4];
IJth(J, 122, 121) =  -p[4];
IJth(J, 122, 122) =  -p[4];
IJth(J, 123, 113) =  p[4];
IJth(J, 123, 123) =  -p[4]-stm*p[5]*p[7];
IJth(J, 123, 124) =  p[8];
IJth(J, 124, 114) =  p[4];
IJth(J, 124, 123) =  stm*p[5]*p[7];
IJth(J, 124, 124) =  -p[4]-p[8]-p[9];
IJth(J, 125, 115) =  p[4];
IJth(J, 125, 124) =  p[9];
IJth(J, 125, 125) =  -p[4]-p[10];
IJth(J, 126, 8) =  -p[4];
IJth(J, 126, 25) =  p[4];
IJth(J, 126, 117) =  p[4];
IJth(J, 126, 126) =  p[4]*-2.0;
IJth(J, 127, 118) =  p[4];
IJth(J, 127, 126) =  p[4];
IJth(J, 127, 127) =  p[4]*-2.0;
IJth(J, 128, 119) =  p[4];
IJth(J, 128, 127) =  p[4];
IJth(J, 128, 128) =  p[4]*-2.0;
IJth(J, 129, 120) =  p[4];
IJth(J, 129, 128) =  p[4];
IJth(J, 129, 129) =  p[4]*-2.0;
IJth(J, 130, 121) =  p[4];
IJth(J, 130, 129) =  p[4];
IJth(J, 130, 130) =  p[4]*-2.0;
IJth(J, 131, 70) =  p[3];
IJth(J, 131, 122) =  p[4];
IJth(J, 131, 130) =  -p[4];
IJth(J, 131, 131) =  -p[4];
IJth(J, 132, 123) =  p[4];
IJth(J, 132, 132) =  -p[4]-stm*p[5]*p[7];
IJth(J, 132, 133) =  p[8];
IJth(J, 133, 124) =  p[4];
IJth(J, 133, 132) =  stm*p[5]*p[7];
IJth(J, 133, 133) =  -p[4]-p[8]-p[9];
IJth(J, 134, 125) =  p[4];
IJth(J, 134, 133) =  p[9];
IJth(J, 134, 134) =  -p[4]-p[10];
IJth(J, 135, 9) =  -p[4];
IJth(J, 135, 26) =  p[4];
IJth(J, 135, 127) =  p[4];
IJth(J, 135, 135) =  p[4]*-2.0;
IJth(J, 136, 128) =  p[4];
IJth(J, 136, 135) =  p[4];
IJth(J, 136, 136) =  p[4]*-2.0;
IJth(J, 137, 129) =  p[4];
IJth(J, 137, 136) =  p[4];
IJth(J, 137, 137) =  p[4]*-2.0;
IJth(J, 138, 130) =  p[4];
IJth(J, 138, 137) =  p[4];
IJth(J, 138, 138) =  p[4]*-2.0;
IJth(J, 139, 71) =  p[3];
IJth(J, 139, 131) =  p[4];
IJth(J, 139, 138) =  -p[4];
IJth(J, 139, 139) =  -p[4];
IJth(J, 140, 132) =  p[4];
IJth(J, 140, 140) =  -p[4]-stm*p[5]*p[7];
IJth(J, 140, 141) =  p[8];
IJth(J, 141, 133) =  p[4];
IJth(J, 141, 140) =  stm*p[5]*p[7];
IJth(J, 141, 141) =  -p[4]-p[8]-p[9];
IJth(J, 142, 134) =  p[4];
IJth(J, 142, 141) =  p[9];
IJth(J, 142, 142) =  -p[4]-p[10];
IJth(J, 143, 10) =  -p[4];
IJth(J, 143, 27) =  p[4];
IJth(J, 143, 136) =  p[4];
IJth(J, 143, 143) =  p[4]*-2.0;
IJth(J, 144, 137) =  p[4];
IJth(J, 144, 143) =  p[4];
IJth(J, 144, 144) =  p[4]*-2.0;
IJth(J, 145, 138) =  p[4];
IJth(J, 145, 144) =  p[4];
IJth(J, 145, 145) =  p[4]*-2.0;
IJth(J, 146, 72) =  p[3];
IJth(J, 146, 139) =  p[4];
IJth(J, 146, 145) =  -p[4];
IJth(J, 146, 146) =  -p[4];
IJth(J, 147, 140) =  p[4];
IJth(J, 147, 147) =  -p[4]-stm*p[5]*p[7];
IJth(J, 147, 148) =  p[8];
IJth(J, 148, 141) =  p[4];
IJth(J, 148, 147) =  stm*p[5]*p[7];
IJth(J, 148, 148) =  -p[4]-p[8]-p[9];
IJth(J, 149, 142) =  p[4];
IJth(J, 149, 148) =  p[9];
IJth(J, 149, 149) =  -p[4]-p[10];
IJth(J, 150, 11) =  -p[4];
IJth(J, 150, 28) =  p[4];
IJth(J, 150, 144) =  p[4];
IJth(J, 150, 150) =  p[4]*-2.0;
IJth(J, 151, 145) =  p[4];
IJth(J, 151, 150) =  p[4];
IJth(J, 151, 151) =  p[4]*-2.0;
IJth(J, 152, 73) =  p[3];
IJth(J, 152, 146) =  p[4];
IJth(J, 152, 151) =  -p[4];
IJth(J, 152, 152) =  -p[4];
IJth(J, 153, 147) =  p[4];
IJth(J, 153, 153) =  -p[4]-stm*p[5]*p[7];
IJth(J, 153, 154) =  p[8];
IJth(J, 154, 148) =  p[4];
IJth(J, 154, 153) =  stm*p[5]*p[7];
IJth(J, 154, 154) =  -p[4]-p[8]-p[9];
IJth(J, 155, 149) =  p[4];
IJth(J, 155, 154) =  p[9];
IJth(J, 155, 155) =  -p[4]-p[10];
IJth(J, 156, 12) =  -p[4];
IJth(J, 156, 29) =  p[4];
IJth(J, 156, 151) =  p[4];
IJth(J, 156, 156) =  p[4]*-2.0;
IJth(J, 157, 74) =  p[3];
IJth(J, 157, 152) =  p[4];
IJth(J, 157, 156) =  -p[4];
IJth(J, 157, 157) =  -p[4];
IJth(J, 158, 153) =  p[4];
IJth(J, 158, 158) =  -p[4]-stm*p[5]*p[7];
IJth(J, 158, 159) =  p[8];
IJth(J, 159, 154) =  p[4];
IJth(J, 159, 158) =  stm*p[5]*p[7];
IJth(J, 159, 159) =  -p[4]-p[8]-p[9];
IJth(J, 160, 155) =  p[4];
IJth(J, 160, 159) =  p[9];
IJth(J, 160, 160) =  -p[4]-p[10];
IJth(J, 161, 13) =  p[4];
IJth(J, 161, 30) =  -p[4];
IJth(J, 161, 75) =  p[3];
IJth(J, 161, 157) =  p[4];
IJth(J, 161, 161) =  -p[4];
IJth(J, 162, 158) =  p[4];
IJth(J, 162, 162) =  -p[4]-stm*p[5]*p[7];
IJth(J, 162, 163) =  p[8];
IJth(J, 163, 159) =  p[4];
IJth(J, 163, 162) =  stm*p[5]*p[7];
IJth(J, 163, 163) =  -p[4]-p[8]-p[9];
IJth(J, 164, 160) =  p[4];
IJth(J, 164, 163) =  p[9];
IJth(J, 164, 164) =  -p[4]-p[10];
IJth(J, 165, 77) =  p[3];
IJth(J, 165, 162) =  -p[4];
IJth(J, 165, 165) =  -stm*p[5]*p[7];
IJth(J, 165, 166) =  p[8];
IJth(J, 166, 78) =  p[3];
IJth(J, 166, 163) =  -p[4];
IJth(J, 166, 165) =  stm*p[5]*p[7];
IJth(J, 166, 166) =  -p[8]-p[9];
IJth(J, 167, 79) =  p[3];
IJth(J, 167, 164) =  -p[4];
IJth(J, 167, 166) =  p[9];
IJth(J, 167, 167) =  -p[10];
IJth(J, 168, 15) =  -stm*p[5]*p[7];
IJth(J, 168, 16) =  -p[8];
IJth(J, 168, 32) =  stm*p[5]*p[7];
IJth(J, 168, 33) =  p[8];
IJth(J, 168, 168) =  -p[8]-p[9]-stm*p[5]*p[7];
IJth(J, 169, 168) =  p[9];
IJth(J, 169, 169) =  -p[10]-stm*p[5]*p[7];
IJth(J, 169, 170) =  p[8];
IJth(J, 170, 16) =  -p[9];
IJth(J, 170, 33) =  p[9];
IJth(J, 170, 169) =  stm*p[5]*p[7];
IJth(J, 170, 170) =  -p[8]-p[9]-p[10];

 
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
