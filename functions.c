/* CS 140
 * Assignment 2 : Matrix Vector Multiplication and the Power Method 
 * Group members : Jonathon Day , Noel Vargas
 * */

/* This file provides the placeholder function definitions, where you will be
 * implementing the assignment algorithms. You will be required to turn in
 * only this file during the submission, where it will be compiled together
 * with our main function and tested. It is therefore required that you keep the
 * function declaration formats unchanged.
 */

#include "powermethod.h"

// Subroutine for generating the input matrix (just one thread's part)
void generatematrix(double * mat, int size)
{
  int myrank,nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //use these to know what rows are getting
  //To Do: make 1111
  //            2222
  //            3333

  int i;
    int fullSize = size * size;
  for (i = 0; i < fullSize; i++ ){
    *(mat + i) = floor(i / size) + 1.0; // every member of matrix is equal to row number
  }
}

// Subroutine to generate a start vector
void generatevec(double * x,int size)
{
  int i;
  for (i = 0; i < size; i++ ){
    *(x + i) = 1; // vector of 1s
  }
}

// Subroutine for the power method, to return the spectral radius
double powerMethod(double * mat, double * x, int size, int iter)
{
  int n = size;
  double *lambda;
  broadcastVector(x, size);
  int iterCount;

 for (iterCount = 0; iterCount < iter; ++iterCount) {
   double *calculatedValues;
  matVec(mat, x,calculatedValues, n,size);
  gatherNewVec(calculatedValues, n, x);
  broadcastVector(x, size);

  double sum;
  sum = norm2(calculatedValues, n);

  updateLambdaVec(lambda,x,sum,n);

  broadcastVector(x, size);

 }

  return *lambda;
}


void updateLambdaVec(double * lambda, double* x,double sum, int n){
    int count;
    for(count = 0; count < n; ++count){
        *(lambda+count) = *(x+count)/sum;
        *(x+count) = *(lambda+ count);
    }
}








void matVec(double *mat, double *vec, double *local_vec, int nrows, int size){

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  free(local_vec);
    local_vec = (double *) calloc(nrows /nprocs, sizeof(double)); //num of processors out of scope????
  //returnMatrix = new [n/ PROC] double ;
  int rowCount, colCount;
    int rowsToParse = nrows/nprocs;
    for(rowCount = 0; rowCount < rowsToParse; ++rowCount){
      *(local_vec + rowCount) = 0;
      for(colCount = 0; colCount < nrows; ++colCount) {
        *(local_vec + rowCount) += *(mat+colCount) *  *(vec+colCount);
      }

    }
}


void broadcastVector(double * vector,int size){

  // Can include better implementation later.
  MPI_Bcast(&vector, size, MPI_INT, 0, MPI_COMM_WORLD);


}

void gatherNewVec(double *vecValues, int n, double * vec){
  // Can optimize later
  double * totalVec;
  int myrank,nprocs;
  int rowSplit = n/nprocs;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MPI_Gather(&vecValues, rowSplit, MPI_INT, &vec, n, MPI_INT, 0, MPI_COMM_WORLD);

}



void receiveSquares(double * square, double *sum){   //maybe need to change return type
  MPI_Reduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  int myrank,nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

}


double norm2(double *x, int size){

    int myrank,nprocs;
    int n = size;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int finishPoint = (myrank + 1) * (n/nprocs);
    double partialSum = 0.0;
    int i;
    for(i = myrank * (n/nprocs); i < finishPoint; ++i){
        double tempSquare = *(x + i);
        partialSum = partialSum + pow(tempSquare, 2);
    }

    double sum = 0.0;
    MPI_Reduce(&partialSum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    return sum;
}