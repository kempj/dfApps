// Static blocked LU Decomposition

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

#define CHECK 1
#define MEASURE 1

void Print_Matrix(double * const v, const int numBlocks, const int size);
void ProcessDiagonalBlock(double * const A, const int L1, const int size);
void ProcessBlockOnRow(double * const A, const double * const D, const int L1, const int L3, const int size);
void ProcessBlockOnColumn(double * const A, const double * const D, const int L1, const int L2, const int size);
void ProcessInnerBlock(double * const A, const double * const remain, const double * const C, const int L1, const int L2, const int L3, const int size);
void InitMatrix2(double *A, const int size);
void InitMatrix3(double *A, const int size);

void stage1(double * const A, const int offset, const int * const sizedim, const int * const start, const int size, const int numBlocks);
void stage2(double * const A, const int offset, const int * const sizedim, const int * const start, const int size, const int numBlocks);
void stage3(double * const A, const int offset, const int * const sizedim, const int * const start, const int size, const int numBlocks);

void LU(double * const A, const int size, const int numBlocks);
void checkResult(double * const A, const double * const A2, const int size);

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
    if( argc > 3 )
    checkResult( A, A2, size);
    free(A);
    free(A2);
    return 0;
}

void LU(double * const A, const int size, const int numBlocks) {
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
    ProcessDiagonalBlock( &A[offset*size+offset], size-offset, size );
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000);
    free(start);
    free(sizedim);
}

void checkResult(double * const A, const double * const A2, const int size) {
    int i, j, k, temp2, temp = 0;
    double *L = (double*)calloc( size*size, sizeof(double) );
    double *U = (double*)calloc( size*size, sizeof(double) );

    if( L==NULL || U==NULL ) {
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

void stage1(double * const A, const int offset, const int * const sizedim, const int * const start, const int size, const int numBlocks)
{
    ProcessDiagonalBlock(&A[offset*size+offset], sizedim[0], size);
}

void stage2(double * const A, const int offset, const int * const sizedim, const int * const start, const int size, const int numBlocks)
{
    int i, L1, L2, L3;
    L1 = sizedim[0]; 
    for(i=1;i<numBlocks;i++){
        L2 = sizedim[i];
        L3 = sizedim[i];
        ProcessBlockOnColumn( &A[(offset+start[i])*size + offset],
                              &A[ offset          *size + offset],
                              L1, L2, size );
    }
    for(i=1;i<numBlocks;i++){
        L2 = sizedim[i];
        L3 = sizedim[i];
        ProcessBlockOnRow( &A[offset*size+(offset+start[i])],
                           &A[offset*size+offset],
                           L1, L3, size );
    }
}

void stage3(double * const A, const int offset, const int * const sizedim, const int * const start, const int size, const int numBlocks)
{
    int i, j, L1, L2, L3;
    L1 = sizedim[0];
    for(i=1; i < numBlocks; i++)
        for(j=1; j < numBlocks; j++){
            L2 = sizedim[i];
            L3 = sizedim[j];
            ProcessInnerBlock( &A[(offset+start[i])*size + (offset+start[j])],
                               &A[ offset          *size + (offset+start[j])], 
                               &A[(offset+start[i])*size + offset], 
                               L1, L2, L3, size );
        }
}

void ProcessDiagonalBlock(double * const A, const int L1, const int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=i+1;j<L1;j++){
            A[j*size+i]/=A[i*size+i];
            for(k=i+1;k<L1;k++)
                A[j*size+k] = A[j*size+k] - A[j*size+i]*A[i*size+k];
        }
}

void ProcessBlockOnColumn(double * const A, const double * const D, const int L1, const int L2, const int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=0;j<L2;j++){
            A[j*size+i]/=D[i*size+i];
            for(k=i+1;k<L1;k++)
                A[j*size+k]+=-A[j*size+i]*D[i*size+k];
        }
}

void ProcessBlockOnRow(double * const A, const double * const D, const int L1, const int L3, const int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=i+1;j<L1;j++)
            for(k=0;k<L3;k++)
                A[j*size+k]+=-D[j*size+i]*A[i*size+k];
}

void ProcessInnerBlock( double * const A, const double * const remain, const double * const C, const int L1, const int L2, const int L3, const int size)
{
    int i,j,k;
    for(i=0;i<L1;i++)
        for(j=0;j<L2;j++)
            for(k=0;k<L3;k++)
                A[j*size+k]+=-C[j*size+i]*remain[i*size+k];
}

void Print_Matrix (double * const v, const int numBlocks, const int size)
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

void InitMatrix2(double *A, const int size)
{
    long long i, j, k;
    for(k = 0; k < size; k++)
        for(i = k; i < size; i++)
            for(j = k; j < size; j++)
                A[i*size + j] += 1;
}

void InitMatrix3(double *A, const int size)
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
#pragma omp for private(i,j)
        for(i = 0; i < size; i++)
            for(j = 0; j < size; j++){
                if(i >= j)
                    L[i*size + j] = i-j+1;
                else
                    L[i*size + j] = 0;
                if(i <= j)
                    U[i*size + j] = j-i+1;
                else
                    U[i*size + j] = 0;
            }
#pragma omp for private(i,j,k)
        for(i = 0; i < size; i++)
            for(j = 0; j < size; j++)
                for(k = 0; k < size; k++)
                    A[i*size + j] += L[i*size + k] * U[k*size + j];
    }
    free(L);
    free(U);
}
