/* CS 140
 * Assignment 2 : Matrix Vector Multiplication and the Power Method 
 * Group members : <Team-member-1> , <Team-member-2>
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
  int i;
  for (i = 0; i < size; i++ ){
    *(mat + i) = -500 + rand() / RAND_MAX / 1000; //random num between -500 and 500
  }
}

// Subroutine to generate a start vector
void generatevec(double * x,int size)
{
  generatematrix(x, size);
}

// Subroutine for the power method, to return the spectral radius
double powerMethod(double * mat, double * x, int size, int iter)
{
  int n = sqrt (size);
  MPI_Bcast(&x, size, MPI_INT, 0, MPI_COMM_WORLD);
double * newVec;
// CALL newVec = multiply_Matrix();

// USE GATHER TO GET ALL THE NEW VECTOR VALUES


  MPI_Gather(newVec, 1, MPI_INT, x, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // gotta figure out the arguments

  // NEW UPDATE OF VECTOR SO WE CAN SQUARE THE SHIT.
  MPI_Bcast(&x, size, MPI_INT, 0, MPI_COMM_WORLD);

double* square;
  square = squareVector();


// USE REDUCE TO GET THE SQUARES BACK IN SUM FORM
//if proc is 1
  // go ahead and get that new vec/ magnitude

// Next is to get the
  return lambda;
}

double* multiply_Matrix(double * mat, double * x, int n){


  double *returnMatrix;
  returnMatrix = new [n/ PROC] double ;
    for(int rowCount = 0; rowCount < n; ++rowCount){
      *returnMatrix[rowCount] = 0;
      for(int colCount = 0; colCount < n; ++colCount) {
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




}

