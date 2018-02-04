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
  for (i = 0; i < size; i++ ){
    *(mat + i) = floor(i / sqrt(size)) + 1.0; // every member of matrix is equal to row number
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
  int n = sqrt (size);
  double lambda;
  broadcastVector(x, size);
  int iterCount;
 for (iterCount = 0; iterCount < iter; ++iterCount) {
   double *calculatedValues;

  calculatedValues = multiply_Matrix(mat, x, n);
  gatherNewVec(calculatedValues, n, x);
  broadcastVector(x, size);

  double sum;
  sum = norm2(calculatedValues, n );
  lambda = calculatedValues/sum;

  // different for every process.
   updateVecValue(x,sum,n);
  broadcastVector(x, size);

 }

  return lambda;
}



void updateVecValue(double *x, double norm, int n){
  int count;
  for (count = 0; count < n; ++count){
    *(x + count) = (*(x + count)) / norm;
  }

}






double*  multiply_Matrix(double * mat, double * x, int n){


  double *returnMatrix;
  returnMatrix = (double *) calloc(n /nprocs, sizeof(double)) //num of processors out of scope????
  //returnMatrix = new [n/ PROC] double ;
  int rowcount, colCount;
    for(rowCount = 0; rowCount < n; ++rowCount){
      *returnMatrix[rowCount] = 0;
      for(colCount = 0; colCount < n; ++colCount) {
        *returnMatrix[rowCount] += *(mat+colCount) * (x+colCount);
      }

    }

  return returnMatrix;
}






double squareVector(double * x, int n){
  double returnValues;
  returnValues = 0;

  for(int count = 0;count < n; count++){
    returnValues[count] += *(x+count) *(x+count);
  }

  return returnValues;
}





void broadcastVector(double * vector,int size){

  // Can include better implementation later.
  MPI_Bcast(&vector, size, MPI_INT, 0, MPI_COMM_WORLD);



}

void gatherNewVec(double * vecValues, int n, double * vec){
  // Can optimize later
  double * totalVec;
  int myrank,nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Gather(&vecValues, n/nprocs, MPI_INT, vec, nprocs, MPI_INT, 0, comm);

}



void receiveSquares(double * square, double *sum){   //maybe need to change return type
  MPI_Reduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  int myrank,nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (nprocs == 1){
    return sum;
  }
  return 0.0;
}

double norm2(double *x, int size){
  int i;
  double sum = 0;
  double tempSquare;
  for(i = 0; i < size; i++){
    tempSquare = *(x + i);
    sum = sum + pow( tempSquare, 2);
  }
  return sqrt(sum);
}