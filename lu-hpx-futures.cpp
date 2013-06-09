// Static blocked LU Decomposition

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

//#include <hpx/hpx_main.hpp>
//#include <hpx/include/lcos.hpp>
//#include <hpx/include/actions.hpp>
////#include <hpx/include/components.hpp>
//#include <hpx/include/iostreams.hpp>
////#include <hpx/include/compression_snappy.hpp>

#include <vector>
using std::vector;
//using hpx::lcos::future;
//using hpx::lcos::wait;
//using hpx::async;

void LU(vector<double> &A, int size, int numBlocks);
void checkResult(vector<double> &A, vector<double> &A2, int size);
void stage1(vector<double> &A, int offset, vector<int> &sizedim, vector<int> &start, int size, int numBlocks);
void stage2(vector<double> &A, int offset, vector<int> &sizedim, vector<int> &start, int size, int numBlocks);
void stage3(vector<double> &A, int offset, vector<int> &sizedim, vector<int> &start, int size, int numBlocks);
void ProcessDiagonalBlock(vector<double> &A, int L1, int offset, int size);
void ProcessBlockOnColumn(vector<double> &A, int p1, int p2, int L1, int L2, int size);
void ProcessBlockOnRow(vector<double> &A, int p1, int p2, int L1, int L2, int size);
void ProcessInnerBlock(vector<double> &A, int p1, int p2, int p3, int L1, int L2, int L3, int size);
void Print_Matrix(vector<double> &v, int numBlocks, int size);
void InitMatrix2(vector<double> &A, int size);
void InitMatrix3(vector<double> &A, int size);

//HPX_PLAIN_ACTION(ProcessBlockOnColumn, ProcessBlockOnColumn_action);
//HPX_PLAIN_ACTION(ProcessBlockOnRow, ProcessBlockOnRow_action);
//HPX_PLAIN_ACTION(ProcessInnerBlock, innerBlock_action);

unsigned long GetTickCount()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + (tv.tv_usec);
}

int main (int argc, char *argv[])
{
    int size = 100, numBlocks = 1;

    if( argc > 1 )
        size = atoi(argv[1]);
    if( argc > 2 )
        numBlocks = atoi(argv[2]);

    printf("size = %d, numBlocks = %d\n", size, numBlocks);

    vector<double> A(size*size, 0);
    InitMatrix3( A, size );
    vector<double> A2(A);
    LU( A, size, numBlocks);
    checkResult( A, A2, size);
    return 0;
}

void LU(vector<double> &  A, int size, int numBlocks)
{
    double t1, t2;
    int i, j, k, remain;
    int itr = 0, offset = 0, Block = 1;
    vector<int> sizedim(numBlocks, 0);
    vector<int> start(numBlocks, 0);

    remain = size;
    t1 = GetTickCount();
    while(size-offset > numBlocks) {
        for(i=0; i < numBlocks; i++) {
            if(i < remain % numBlocks) {
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
        remain -= sizedim[0];
    }
    //ProcessDiagonalBlock( &A[offset*size+offset], size-offset, size );
    ProcessDiagonalBlock(A, size-offset, offset, size );
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000);
}

void stage1(vector<double> & A,
            int offset, 
            vector<int> & sizedim, 
            vector<int> & start, 
            int size, 
            int numBlocks)
{
    //ProcessDiagonalBlock(&A[offset*size+offset], sizedim[0], size);
    ProcessDiagonalBlock(A, sizedim[0], offset, size);
}

void stage2(vector<double> & A, 
            int offset, 
            vector<int> & sizedim, 
            vector<int> & start, 
            int size, 
            int numBlocks)
{
    int L1, L2, L3;
    L1 = sizedim[0]; 
    //ProcessBlockOnColumn_action blockCol;
    //ProcessBlockOnRow_action blockRow;
    //future<void>* futures = new future<void>[2*numBlocks];
    for(int i=1;i<numBlocks;i++){
        L2 = sizedim[i];
        L3 = sizedim[i];
        //futures[2*i] = hpx::async(blockCol, hpx::find_here(),
        //                              &A[(offset+start[i])*size + offset],
        //                              &A[ offset          *size + offset],
        //                              L1, L2, size ) ;
        //ProcessBlockOnColumn( &A[(offset+start[i])*size + offset],
        //                      &A[ offset          *size + offset],
        //                      L1, L2, size );
        ProcessBlockOnColumn( A, 
                              (offset+start[i])*size + offset, 
                              offset*size + offset,
                              L1, L2, size );
        //futures[2*i+1] = hpx::async(blockRow, hpx::find_here(),
        //                              &A[offset*size+(offset+start[i])],
        //                              &A[offset*size+offset],
        //                              L1, L3, size ) ;
        //ProcessBlockOnRow( &A[offset*size+(offset+start[i])],
        //                   &A[offset*size+offset],
        //                   L1, L3, size );
        ProcessBlockOnRow( A,
                           offset*size+(offset+start[i]),
                           offset*size+offset,
                           L1, L3, size );
    }
    //vector<future<void>> f(&futures[0], &futures[0] + numBlocks);
    //wait(f);
}

void stage3(vector<double> &A, 
            int offset, 
            vector<int> &sizedim,
            vector<int> &start, 
            int size, 
            int numBlocks)
{
    int L1, L2, L3;
    L1 = sizedim[0];

    for(int i=1; i < numBlocks; i++)
        for(int j=1; j < numBlocks; j++){
            L2 = sizedim[i];
            L3 = sizedim[j];
            ProcessInnerBlock(A, (offset+start[i])*size + (offset+start[j]),
                                  offset          *size + (offset+start[j]), 
                                 (offset+start[i])*size + offset, 
                               L1, L2, L3, size );
        }
}

void ProcessDiagonalBlock(vector<double> &A, int L1, int offset, int size)
{
    int pivot = offset*size + offset;
    for(int i = 0; i < L1; i++)
        for(int j = i+1; j<L1; j++){
            A[pivot+j*size+i] /= A[pivot+i*size+i];
            for(int k = i+1; k < L1; k++)
                A[pivot+j*size+k] = A[pivot+j*size+k] - A[pivot+j*size+i] * A[pivot+i*size+k];
        }
}

void ProcessBlockOnColumn(vector<double> &A,
                          int p1, int p2,
                          int L1, int L2, int size)
{
    for(int i=0; i < L1; i++)
        for(int j=0; j < L2; j++){
            A[p1+j*size+i] /= A[p2+i*size+i];
            for(int k = i+1; k < L1; k++)
                A[p1+j*size+k] += -A[p1+j*size+i] * A[p2+i*size+k];
        }
}

void ProcessBlockOnRow(vector<double> &A, 
                       int p1, int p2, 
                       int L1, int L3, int size)
{
    for(int i=0;i<L1;i++)
        for(int j=i+1;j<L1;j++)
            for(int k=0;k<L3;k++)
                A[p1 + j*size+k]+=-A[p2 + j*size+i]*A[p1+i*size+k];
}

void ProcessInnerBlock(vector<double> &A, 
                       int p1, int p2, int p3,
                       int L1, int L2, int L3, int size)
{
    for(int i=0; i < L1; i++)
        for(int j=0; j < L2; j++)
            for(int k=0; k < L3; k++)
                A[p1+j*size+k] += -A[p3+j*size+i] * A[p2+i*size+k];
}

void checkResult(vector<double> &A, vector<double> &A2, int size) 
{
    int i, j, k, temp2, temp = 0;
    vector<double> L(size*size, 0);
    vector<double> U(size*size, 0);

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
            if((A2[i*size+j]-temp2)/A2[i*size+j] > 0.1 || (A2[i*size+j]-temp2)/A2[i*size+j] < -0.1 )
                temp++;
        }
    printf("Errors = %d \n", temp);
}

void Print_Matrix(vector<double> &v, int numBlocks, int size)
{
    printf( "\n" );
    for(int i = 0; i < numBlocks; i++){
        for(int j = 0; j < size; j++)
            printf( "%.2f,", v[i*size + j] );
        printf( "\n" );
    }
    printf( "\n" );
}

void InitMatrix2(vector<double> &A, int size)
{
    for(long long k = 0; k < size; k++)
        for(long long i = k; i < size; i++)
            for(long long j = k; j < size; j++)
                A[i*size + j] += 1;
}

void InitMatrix3(vector<double> &A, int size)
{
    vector<double> L(size*size, 0), U(size*size,0);
    //TODO: parallelize?
    for(long long i = 0; i < size; i++)
        for(long long j = 0; j < size; j++){
            if(i >= j)
                L[i*size + j] = i-j+1;
            else
                L[i*size + j] = 0;
            if(i <= j)
                U[i*size + j] = j-i+1;
            else
                U[i*size + j] = 0;
        }
    for(long long i = 0; i < size; i++)
        for(long long j = 0; j < size; j++)
            for(long long k = 0; k < size; k++)
                A[i*size + j] += L[i*size + k] * U[k*size + j];
}
