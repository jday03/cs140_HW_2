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
    int amountOfColumns;
    double *tempMat;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int i, j;
    int x;
    amountOfColumns = size / nprocs;

    tempMat = (double *) calloc(size * size, sizeof(double));
    for(i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (j < i + 1) {
                x = i * size + j;
                tempMat[x] = i + 1;
            } else {
                x = i * size + j;
                tempMat[x] = 0;
            }

        }
    }

    for(i=0; i< size*size; i++){
       // printf("tempMat[%d] = %f\n" , i, tempMat[i]);
    }
    for(i=0; i < size * amountOfColumns; i++){
        mat[i] = tempMat[i + myrank*amountOfColumns * size] ;
    }

  //  printf("myrank = %d\n", myrank);
    for(i = 0; i < size *amountOfColumns; i++){
       // printf("mat[%d] : %f\n", i, mat[i]);

    }

    free(tempMat);
}


// Subroutine to generate a start vector
void generatevec(double * x,int size)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank == 0) {
        int i = 0;
        double *iter;
        for (iter = x; i < size; iter++) {
            x[i] = 1.0;
           // printf("VEC %f\n", x[i]);
            i++;
        }
    }
}

// Subroutine for the power method, to return the spectral radius
double powerMethod(double * mat, double * x, int size, int iter)
{
    int myrank,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double *lambda = (double *) calloc(1, sizeof(double));
    double *calculatedValues= (double *) calloc(size/nprocs, sizeof(double));




  broadcastVector(x, size);
    int counter;

  int iterCount;

    for (iterCount = 0; iterCount < iter; ++iterCount) {
  matVec(mat, x,calculatedValues, (size /nprocs),size);


     MPI_Allgather(calculatedValues,2*size/nprocs,MPI_INT, x, 2*size/nprocs,MPI_INT,MPI_COMM_WORLD);




                    int counter;
            for(counter = 0; counter < size;++counter) {
               // printf("FULL VEC:item %u : %f \n", counter, *(x+ counter));
            }

       // printf("PROCSUMMMMM333 %u vector Num 7: %f \n", myrank, (*(x + 7)));


   broadcastVector(x, size);

   double sum;
   sum = norm2(x, size);

        if(myrank == 0){
        //printf("SUMMATION AT END: %f",sum);
        sum = sqrt(sum);
            updateLambdaVec(lambda,x,sum,size);
            *lambda = sum;
        }


   broadcastVector(x, size);

 }

  return *lambda;
}


void updateLambdaVec(double * lambda, double* x,double sum, int n){
    int count;
    for(count = 0; count < n; ++count){
        *(x+count) = *(x+count)/sum;
    }
}








void matVec(double *mat, double *vec, double *local_vec, int nrows, int size){

    int nprocs;
    int myrank = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


  int rowCount, colCount;
    int rowsToParse = nrows;
    for(rowCount = 0; rowCount < rowsToParse; ++rowCount){
        *(local_vec + rowCount) = 0;
        double sum = 0;

        for(colCount = 0; colCount < size; ++colCount) {
          //  printf("ROWCOUNT IS %u\n", rowCount);

            *(local_vec + rowCount) += *(mat+((size * rowCount)+colCount)) *  *(vec+colCount);

            sum += *(mat+colCount) *  *(vec+colCount);

        }
       // printf("PROC %u vector Num %u: %f should be %f \n", myrank,rowCount, *(local_vec + rowCount), sum);

    }
}


void broadcastVector(double * vector,int size){
//printf("SIZE %d", sizeof(double));
  // Can include better implementation later.
  MPI_Bcast(vector, 2* size, MPI_INT, 0, MPI_COMM_WORLD);


}

void gatherNewVec(double *vecValues, int n, double * vec){
  // Can optimize later
  double *totalVec;
  int myrank,nprocs;
  int rowSplit = n/nprocs;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    if(myrank == 0){

        MPI_Gather(&vecValues, n/nprocs, MPI_INT, &vec, n/nprocs, MPI_INT, 0, MPI_COMM_WORLD);


    }
    else MPI_Gather(&vecValues, n/nprocs, MPI_INT, NULL, n/nprocs, MPI_INT, 0, MPI_COMM_WORLD);

}



void receiveSquares(double * square, double *sum){   //maybe need to change return type
  MPI_Reduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  int myrank,nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

}


double norm2(double* x, int size){


    int myrank,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    int finishPoint = (myrank + 1) * (size/nprocs);
    double partialSum = 0.0;
    int i;
    for(i = myrank * (size/nprocs); i < finishPoint; ++i){
        double tempSquare = *(x + i);
        partialSum = partialSum + pow(tempSquare, 2);

    }
   // printf("PROCSUMMMMM %u vector partialSum = %f \n", myrank, partialSum);

    double sum =0.0;
    MPI_Reduce(&partialSum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //printf("SUM after reduce: %u vector Sum = %f \n", myrank, sum);

    return sum;
}