// Static blocked LU Decomposition

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

//#define DEBUG 1
#define CHECK 1
#define MEASURE 1

void Print_Matrix(double *v, int numBlocks, int size);
void ProcessDiagonalBlock(double *A, int L1, int size);
void ProcessBlockOnRow(double *A, double *D, int L1, int L3, int size);
void ProcessBlockOnColumn(double *A, double *D, int L1, int L2, int size);
void ProcessInnerBlock(double *A, double *remain, double *C, int L1, int L2, int L3, int size);
void InitMatrix2(double *A, int size);
void InitMatrix3(double *A, int size);

void stage1(double *A, int offset, int *sizedim, int *start, int size, int numBlocks);
void stage2(double *A, int offset, int *sizedim, int *start, int size, int numBlocks);
void stage3(double *A, int offset, int *sizedim, int *start, int size, int numBlocks);

void LU(double *A, int size, int numBlocks);
void checkResult(double *A, double *A2, int size);

unsigned long GetTickCount()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + (tv.tv_usec);
}

int main (int argc, char *argv[])
{
    double *A, *A2;
    int size = 100, numBlocks = 1;

    if( argc > 1 )
        size = atoi(argv[1]);
    if( argc > 2 )
        numBlocks = atoi(argv[2]);

    printf("size = %d, numBlocks = %d\n", size, numBlocks);
    A  = (double*)calloc( size*size, sizeof(double) );
    A2 = (double*)calloc( size*size, sizeof(double) );

    if(A==NULL || A2==NULL) {
        printf( "Can't allocate memory\n" );
        exit(1);
    }

    InitMatrix3( A, size );

    memcpy( A2, A, size*size*sizeof(double) );
    LU( A, size, numBlocks);
#ifdef CHECK
    checkResult( A, A2, size);
#endif
    free(A);
    free(A2);
    return 0;
}

void LU(double *A, int size, int numBlocks) {
    double t1, t2;
    int i, j, k, remain;
    int itr = 0, offset = 0, Block = 1;
    int *sizedim = (int*)malloc( numBlocks*sizeof(int) );
    int *start   = (int*)malloc( numBlocks*sizeof(int) );
    if(start == NULL || sizedim == NULL ) {
        printf( "Can't allocate memory\n" );
        exit(1);
    }

    remain = size;
    t1 = GetTickCount();
#pragma omp parallel
    {
#pragma omp master
        {
            while(size-offset > numBlocks){
                for(i=0; i < numBlocks; i++){
                    if(i < remain%numBlocks){
                        sizedim[i] = remain/numBlocks+1;
                        start[i]   = (remain/numBlocks+1)*i;
                    } else {
                        sizedim[i] = remain/numBlocks;
                        start[i]   = (remain/numBlocks+1)*(remain%numBlocks) +
                                     (remain/numBlocks)*(i-remain%numBlocks);
                    }
                }
                stage1( A, offset, sizedim, start, size, numBlocks );
                stage2( A, offset, sizedim, start, size, numBlocks );
                stage3( A, offset, sizedim, start, size, numBlocks );
                offset += sizedim[0];
                remain = remain-sizedim[0];
            }
        }
    }
    ProcessDiagonalBlock( &A[offset*size+offset], size-offset, size );
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000);
    free(start);
    free(sizedim);
}

void checkResult(double *A, double *A2, int size) {
    int i, j, k, temp2, temp = 0;
    double *L = (double*)calloc( size*size, sizeof(double) );
    double *U = (double*)calloc( size*size, sizeof(double) );

    if(L== NULL || U==NULL ) {
        printf( "Can't allocate memory\n" );
        exit(1);
    }

    for(i=0;i<size;i++)
        for(j=0;j<size;j++)
            if (i>j)
                L[i*size+j] = A[i*size+j];
            else
                U[i*size+j] = A[i*size+j];
    for(i=0;i<size;i++)
        L[i*size+i] = 1;

    for(i=0;i<size;i++)
        for(j=0;j<size;j++){
            temp2=0;
            for(k=0;k<size;k++)
                temp2+=L[i*size+k]*U[k*size+j];
            if( (A2[i*size+j]-temp2)/A2[i*size+j] > 0.1 || 
                (A2[i*size+j]-temp2)/A2[i*size+j] < -0.1 )
                temp++;
        }
    printf("Errors = %d \n", temp);
    free(L);
    free(U);
}

// always start[0] is 0, so it is not used;
void stage1(double *A, int offset, int *sizedim, int *start, int size, int numBlocks)
{
    ProcessDiagonalBlock(&A[offset*size+offset], sizedim[0], size);
}

void stage2(double *A, int offset, int *sizedim, int *start, int size, int numBlocks)
{
    int B, i, L1, L2, L3;
    B = L1 = sizedim[0]; 
    /* Processing only one big block in column and row */
    for(i=1;i<numBlocks;i++){
        L2 = sizedim[i];
        L3 = sizedim[i];
#pragma omp task firstprivate(i, L1, L2, offset, size)
        ProcessBlockOnColumn( &A[(offset+start[i])*size + offset],
                              &A[ offset          *size + offset],
                              L1, L2, size );
#pragma omp task firstprivate(i, L1, L2, offset, size)
        ProcessBlockOnRow( &A[offset*size+(offset+start[i])],
                           &A[offset*size+offset],
                           L1, L3, size );
    }
#pragma omp taskwait
}

void stage3(double *A, int offset, int *sizedim, int *start, int size, int numBlocks)
{
    int i, j, B, L1, L2, L3;
    B = L1 = sizedim[0];
    for(i=1; i < numBlocks; i++)
        for(j=1; j < numBlocks; j++){
            L2 = sizedim[i];
            L3 = sizedim[j];
#pragma omp task firstprivate(i,j,numBlocks,size,offset,L1,L2,L3)
            ProcessInnerBlock( &A[(offset+start[i])*size + (offset+start[j])],
                               &A[ offset          *size + (offset+start[j])], 
                               &A[(offset+start[i])*size + offset], 
                               L1, L2, L3, size );
        }
#pragma omp taskwait
}

/* *A is a pointer to the block processed 
 * The size of the diagonal block is L1xL1 
 * size is the size of the matrix in one dimension 
 */
void ProcessDiagonalBlock(double *A, int L1, int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=i+1;j<L1;j++){
            A[j*size+i]/=A[i*size+i];
            // DAXPY(&A[j*size+(i+1)],&A[i*size+(i+1)] ,-A[j*size+i],L1-(i+1),1);
            for(k=i+1;k<L1;k++)
                A[j*size+k] = A[j*size+k] - A[j*size+i]*A[i*size+k];
        }
}

/* *A is a pointer to the column block processed 
 * *D is a pointer to the diagonal block required 
 * The size of the column block is L2xL1 
 * The size of the diagonal block is L1xL1 
 */
void ProcessBlockOnColumn(double *A, double *D, int L1, int L2, int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=0;j<L2;j++){
            A[j*size+i]/=D[i*size+i];
            //  DAXPY(&A[j*size+(i+1)],&D[i*size+(i+1)],-A[j*size+i],L1-(i+1),1);
            for(k=i+1;k<L1;k++)
                A[j*size+k]+=-A[j*size+i]*D[i*size+k];
        }
}

/* *A is a pointer to the row block processed 
 * *D is a pointer to the diagonal block required 
 * The size of the row block is L1xL3 
 * The size of the diagonal block is L1xL1 
 */
void ProcessBlockOnRow(double *A, double *D, int L1, int L3, int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=i+1;j<L1;j++)
            // DAXPY(&A[size*j],&A[size*i],-D[j*size+i],L3,1);
            for(k=0;k<L3;k++)
                A[j*size+k]+=-D[j*size+i]*A[i*size+k];
}

/* *A is a pointer to the inner block processed 
 * *remain is a pointer to the row block required 
 * *C is a pointer to the column block required 
 * The size of the row block is L1xL3 
 * The size of the column block is L2xL1 
 * The size of the inner block is L2xL3 
 */
void ProcessInnerBlock(double *A, double *remain, double *C, int L1, int L2, int L3, int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=0;j<L2;j++)
            // DAXPY(&A[size*j],&remain[size*i],-C[j*size+i],L3,1);
            for(k=0;k<L3;k++)
                A[j*size+k]+=-C[j*size+i]*remain[i*size+k];
}

void Print_Matrix (double *v, int numBlocks, int size)
{
    int i,j;
    printf( "\n" );
    for(i = 0; i < numBlocks; i++){
        for(j = 0; j < size; j++)
            printf( "%.2f,", v[i*size + j] );
        printf( "\n" );
    }
    printf( "\n" );
}

void InitMatrix2(double *A, int size)
{
    long long i, j, k;
    //for(i = 0; i < size*size; i++) A[i] = 0;
    for(k = 0; k < size; k++)
        for(i = k; i < size; i++)
            for(j = k; j < size; j++)
                A[i*size + j] += 1;
}

void InitMatrix3(double *A, int size)
{
    long long i, j, k;
    double *L, *U;
    L = (double*)calloc(size*size, sizeof(double));
    U = (double*)calloc(size*size, sizeof(double));
    if(L == NULL || U == NULL ) {
        printf( "Can't allocate memory\n" );
        exit(1);
    }
#pragma omp parallel 
    {
#pragma omp forprivate(i,j)
        for(i = 0; i < size; i++)
            for(j = 0; j < size; j++){
                //A[i*size + j] = 0;
                if(i >= j)
                    L[i*size + j] = i-j+1;
                else
                    L[i*size + j] = 0;
                if(i <= j)
                    U[i*size + j] = j-i+1;
                else
                    U[i*size + j] = 0;
            }
#pragma omp forprivate(i,j,k)
        for(i = 0; i < size; i++)
            for(j = 0; j < size; j++)
                for(k = 0; k < size; k++)
                    A[i*size + j] += L[i*size + k] * U[k*size + j];
    }
    free(L);
    free(U);
}

    /*   /\* LU DECOMPOSITION *\/ */
    /*   for(k=0;k<size-1;k++){ */
    /*     for(i=k+1;i<size;i++){ */
    /*       A[i*size+k] = A[i*size+k]/A[k*size+k]; */
    /*    /\*  for(i=k+1;i<size;i++) *\/ */
    /*       for(j=k+1;j<size;j++) */
    /*         A[i*size+j] = A[i*size+j] - A[i*size+k]*A[k*size+j]; */
    /*     } */
    /*   } */

/*void InitMatrix(double *A, int size)
  {
  long long	i, j;
  struct timeval	InitTime;

  gettimeofday(&InitTime, NULL);

  srand(InitTime.tv_sec * 1000000 + InitTime.tv_usec);

  for(i = 0; i < size; i++) {
  for(j = 0; j < size; j++) {
  A[i * size + j] = ( rand() - remainAsizeD_numBlocksAX/2) / (remainAsizeD_numBlocksAX/1000 );
  if (i == j) {
  A[i * size + i] *= 10;
  }
//A[i*size+j]=i+j+1;
}
}
} 
*/
