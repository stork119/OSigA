/*
 * -----------------------------------------------------------------
 *  * -----------------------------------------------------------------
 */

/* libraries */
#include "model_rcpp_libraries.h"
#include "model_rcpp_models.h"


/* Problem Constants */

#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL RCONST(1.0e-4) 

/* Functions Called by the Solver */



struct solverData {
    std::vector< std::vector<double> > trajectory;  
    int solver_flag;
};


/*static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);


static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
*/



/* Private functions to output results */
// static void writeToFile(N_Vector y, FILE* outputFile);

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
    int neq;
    int npar;
    double t0;
    double t1;
    double tmult;
    int tnout;
    double t_stm_start;
    double t_stm_end;
};

/*
 *---------------------------------
 * Run odesolver
 *---------------------------------
 */
solverData run_solver(
  double* parameters, double* variables, double stm, int neq, int npar, double t0, double t1, double tmult, int tnout, double t_stm_start, double t_stm_end)
{

  static int (*f)(realtype t, N_Vector yith, N_Vector ydot, void *_user_data);
  static int (*Jac)(long int N, realtype t,
               N_Vector yvec, N_Vector fy, DlsMat J, void *_user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  Jac = Jac_jakstat_means;
  f = f_jakstat_means; 


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
  y = N_VNew_Serial(neq);
  if (check_flag((void *)y, "N_VNew_Serial", 0)){
    return(solver_results);
  }
  abstol = N_VNew_Serial(neq); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)){
    return(solver_results);
  }
  
  std::vector<double> yvector;  
  for(int i = 0; i < neq; ++i){
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
   * user's right hand side function in y'=f(t,y), the inital time t0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, t0, y);
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
  flag = CVDense(cvode_mem, neq);
  if (check_flag(&flag, "CVDense", 1)){
    return(solver_results);
  }

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)){
    return(solver_results);
  }

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when tnout preset output times have been reached.  */
  solver_results.solver_flag = 1;
  iout = 0;
  tout = t1;
  while(1) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    std::vector<double> yvector;  
    for(int i = 0; i < neq; ++i){
      yvector.push_back(Ith(y, i+1));
//      printf("%lf\t", Ith(y, i+1));
    } 
    results.push_back(yvector);
//    printf("\n");

    if (check_flag(&flag, "CVode", 1)){
      solver_results.solver_flag = 0;
      break;
    }
    if (flag == CV_SUCCESS) {
      iout++;
      tout = tout + tmult;
    }
    if (iout == tnout) break;
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
    thread_args->stm,
    thread_args->neq,
    thread_args->npar,
    thread_args->t0,
    thread_args->t1,
    thread_args->tmult,
    thread_args->tnout,
    thread_args->t_stm_start,
    thread_args->t_stm_end 
  );

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
  int model_type,
  double stm,
  double t0,
  double t1,
  double tmult,
  int tnout,
  double t_stm_start,
  double t_stm_end,
  unsigned long time_interval,
  unsigned long time_computation
)
{
  struct modelData model_data = returnModelStructure_jakstat_means();

  //clock_t start = clock();
  Rcpp::List list;
  list["success"] = 0;

  if( parameters.size() != model_data.npar ){
    return(list);

  }
  double parametersArray[model_data.npar + 1];
  for(int i = 0; i < model_data.npar; ++i){
    parametersArray[i+1] = parameters[i];
  } 

  if( variables.size() < model_data.neq ){
    return(list);
  }
  double variablesArray[model_data.neq + 1];
  for(int i = 0; i < model_data.neq; ++i){
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
    thread_args.neq = model_data.neq;
    thread_args.npar = model_data.npar;
    thread_args.t0 = t0;
    thread_args.t1 = t1;
    thread_args.tmult = tmult;
    thread_args.tnout = tnout;   
    thread_args.t_stm_start = t_stm_start;    
    thread_args.t_stm_end = t_stm_end;    


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
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */


/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

// static void writeToFile(N_Vector y, FILE* outputFile)
// {
//   for(int i = 1; i < NEQ; ++i)
//   {
//     fprintf(outputFile, "%lf,", Ith(y, i) );
//   }
//   fprintf(outputFile, "%lf\n", Ith(y, NEQ));
// }

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





