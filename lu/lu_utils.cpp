#include "lu_utils.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <future>

using std::future;
using std::vector;
using std::async;


void initA( vector<double> &L, vector<double> &U, int i, int size)
{
    for(int j = 0; j < size; j++)
        for(int k = 0; k < size; k++)
            A[i * size + j] += L[i * size + k] * U[k*size + j];
}

void initLU( vector<double> &L, vector<double> &U, int i, int size)
{
    for(int j = 0; j < size; j++){
        if(i >= j)
            L[i*size + j] = i-j+1;
        else
            L[i*size + j] = 0;
        if(i <= j)
            U[i*size + j] = j-i+1;
        else
            U[i*size + j] = 0;
    }
}

void InitMatrix3( int size )
{
    vector<double> L, U;
    L.reserve(size*size);
    U.reserve(size*size);
    vector<future<void>> futures, LUfutures;
    futures.reserve(size);
    LUfutures.reserve(size);
    for(int i = 0; i < size; i++) {
        LUfutures.push_back( async( std::launch::async, &initLU, std::ref(L), std::ref(U), i, size));
    }
    for(int i = 0; i < LUfutures.size(); i++) {
        LUfutures[i].wait();
    }
    for(int i = 0; i < size; i++) {
        futures.push_back( async( std::launch::async, &initA, std::ref(L), std::ref(U), i, size));
    }
    for(int i = 0; i < futures.size(); i++) {
        futures[i].wait();
    }
//    hpx::lcos::wait(futures);
}

void Print_Matrix(vector<double> &v, int size)
{
    printf( "\n" );
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++)
            printf( "%5.1f, ", v[i*size + j] );
        printf( "\n" );
    }
    printf( "\n" );
}

void checkResult( vector<double> &originalA, int size )
{
    int errors = 0;
    double temp2;
    vector<double> L(size*size, 0);
    vector<double> U(size*size, 0);
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            if (i>j)
                L[i*size+j] = A[i*size+j];
            else
                U[i*size+j] = A[i*size+j];
    for(int i=0;i<size;i++)
        L[i*size+i] = 1;

    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++){
            temp2=0;
            for(int k=0;k<size;k++)
                temp2+=L[i*size+k]*U[k*size+j];
            if( (originalA[i*size+j]-temp2) / originalA[i*size+j] > 0.1 || (originalA[i*size+j]-temp2) / originalA[i*size+j] < -0.1 ){
                printf("error:[%d][%d] ", i, j);
                errors++;
            }
        }
    if(errors > 0){
        printf("A:\n");
        Print_Matrix(A, size);
        printf("originalA:\n");
        Print_Matrix(originalA, size);
    }

    printf("Errors = %d \n", errors);
}

unsigned long GetTickCount()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + (tv.tv_usec);
}
