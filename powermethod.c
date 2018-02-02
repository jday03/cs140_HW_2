/* CS 140
 * Assignment 2 : Matrix Vector Multiplication and the Power Method
 * */

/* This is a sample main function for the assignment provided for testing purposes.
 * You may make any changes you want to this file for testing purposes. 
 * The function definitions are to be written in functions.c.
 */

#include "powermethod.h"

int main(int argc, char **argv)
{
  double *mat, *vec;
  double spectral_radius;
  double start,stop;

  MPI_Init(&argc,&argv);

  // Get the size of the matrix and number of iterations from the command line arguments.

  int size = atoi(argv[1]);
  int iter = atoi(argv[2]);

  // Get the number of threads and the rank of the current thread.
  
  int myrank,nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // Allocate memory for the matrix (just the part on this thread) and the vector.
  
  mat = (double *) calloc(((size+nprocs-1)/nprocs) * size, sizeof(double));
  vec = (double *) calloc(size, sizeof(double));

  // Thread 0 does all the output.
  
  if( myrank == 0) {
    printf("\nThis program generates a matrix of dimensions %d x %d and runs %d iterations\n", size, size, iter);
    printf("of the power method on it.  For the example matrix in the problem, the answer is %d.\n", size);
    printf("size of mat alloc = %d\n", ( ((size+nprocs-1) / nprocs)*size ) );



  }

  generatematrix(mat, size); // every thread generates its part of the matrix.
  generatevec(vec, size);    // every thread generates all of the start vector.

  // Power method to generate the spectral radius of the matrix.
  
  start = MPI_Wtime();
  spectral_radius = powerMethod(mat, vec, size, iter);
  stop = MPI_Wtime();

  // Print the radius and execution time using thread 0.
  // You can change the print format anyway you want.

  if(myrank == 0) {
    printf("Spectral radius: %lf\n",spectral_radius);
    printf("Time in seconds: %lf\n", stop - start);
  }





  free(mat);
  free(vec);

  MPI_Finalize();
}
