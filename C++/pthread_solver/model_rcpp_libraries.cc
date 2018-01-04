#include "model_rcpp_libraries.h"
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





