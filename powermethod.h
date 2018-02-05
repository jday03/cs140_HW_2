#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

// Function declarations

void generatematrix(double * mat, int size);
void generatevec(double * x, int size);
double powerMethod(double * mat, double * vec, int size, int iter);
double norm2(double *x, int size);
void matVec(double *mat, double *vec, double *local_vec, int nrows, int size);
void updateLambdaVec(double * lambda, double* x,double sum, int n);
double squareVector(double * x, int n);
void broadcastVector(double * vector,int size);
void gatherNewVec(double *vecValues, int n, double * vec);
void receiveSquares(double * square, double *sum);



