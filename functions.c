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
  return lambda;
}
