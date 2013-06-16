// Static blocked LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <hpx/hpx_main.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/actions.hpp>

#include <vector>

using std::vector;
using hpx::lcos::future;
using hpx::lcos::wait;
using hpx::async;

struct block {
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & size;
            ar & start;
            ar & height;
        }
    int size;
    int start;
    int height;
    block(int size, int startAddress, int H) : size(size), start(startAddress), height(H){}
    block() : size(0), start(0), height(0){}
};
void LU( int size, int numBlocks);
void checkResult( vector<double> &A2, int size);

void ProcessDiagonalBlock( int size, block B);
void ProcessBlockOnColumn( int size, block B1, block B2);
void ProcessBlockOnRow( int size, block B1, block B2);
void ProcessInnerBlock( int size, block B1, block B2, block B3);

void getBlockList(vector<vector<block>> &blocks, int size, int numBlocks);

void Print_Matrix(vector<double> &v, int numBlocks, int size);
void InitMatrix3( int size);
void initLoop(int i, int size);

HPX_PLAIN_ACTION( ProcessBlockOnColumn, column_action);
HPX_PLAIN_ACTION( ProcessBlockOnRow, row_action);
HPX_PLAIN_ACTION( ProcessInnerBlock, innerBlock_action);
HPX_PLAIN_ACTION( initLoop, init_action);
HPX_PLAIN_ACTION( ProcessDiagonalBlock, diag_action);

vector<double> A;
vector<double> L;
vector<double> U;
unsigned long GetTickCount()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + (tv.tv_usec);
}

int main (int argc, char *argv[])
{
    unsigned long t1, t2;
    int size = 100, numBlocks = 1;
    if( argc > 1 )
        size = atoi(argv[1]);
    if( argc > 2 )
        numBlocks = atoi(argv[2]);
    printf("size = %d, numBlocks = %d\n", size, numBlocks);
    A.resize(size*size, 0);
    L.resize(size*size, 0);
    U.resize(size*size, 0);
    t1 = GetTickCount();
    InitMatrix3( size );
    t2 = GetTickCount();
    vector<double> A2;
    A2.reserve(size*size);
    for(int i = 0; i < size * size; i++)
        A2[i] = A[i];
    printf("init done, time = %f\n", (t2-t1)/1000000.0);
    t1 = GetTickCount();
    if(numBlocks > 1)
        ProcessDiagonalBlock(size, block(size, 0, size));
    else
        LU( size, numBlocks);
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000.0);
    checkResult( A2, size);
    return 0;
}

void LU( int size, int numBlocks)
{
    vector<vector<block>> blockList;
    column_action blockCol;
    row_action blockRow;
    innerBlock_action blockInner;
    diag_action blockDiag;

        getBlockList(blockList, size, numBlocks);

        for(int i = 0; i < numBlocks; i++) {
            ProcessDiagonalBlock( size, blockList[i][i] );

            for(int j=i + 1; i < numBlocks; j++){
                dataflow<blockCol>( hpx::find_here(), size, blocks[j][i], blocks[i][i] );
                dataflow<blockRow>( hpx::find_here(), size, blocks[i][j], blocks[i][i] );
            }

            for(int j=i+1; i < numBlocks; j++) {
                for(int k=i+1; k < numBlocks; k++) {
                    dataflow<blockInner>( hpx::find_here(), size, blocks[j][k], blocks[i][k], blocks[j][i] );
                }
            }
        }
    } else {
    }
}

void getBlockList(vector<vector<block>> &blockList, int size, int numBlocks)
{
    int blockSize, start, height;
    for(int i=0; i < numBlocks; i++) 
        blockList.push_back(vector<block>());

    height = size/numBlocks;
    if(size%numBlocks > 0)
        height += 1;
    for(int i=0; i < numBlocks; i++) {
        if(i < size % numBlocks) {
            blockSize = size/numBlocks+1;
            start = (size/numBlocks+1)*i;
        } else {
            blockSize = size/numBlocks;
            start = (size/numBlocks+1)*(size%numBlocks) + (size/numBlocks)*(i-size%numBlocks);
        }
        blockList[0].push_back( block( blockSize, start, height));
    }
    for(int i = 1; i < numBlocks; i++) {
        height = blockList[0][i].size;
        for(int j = 0; j < numBlocks; j++) {
            blockSize = blockList[0][j].size;
            start = blockList[i-1][j].start + blockList[i-1][0].height * size;
            blockList[i].push_back( block( blockSize, start, height));
        }
    }
}

void ProcessDiagonalBlock( int size, block B)
{
    for(int i = 0; i < B.size; i++) {
        for(int j = i+1; j < B.size; j++){
            A[B.start+j*size+i] /= A[B.start+i*size+i];
            for(int k = i+1; k < B.size; k++) {
                A[B.start+j*size+k] -= A[B.start+j*size+i] * A[B.start+i*size+k];
            }
        }
    }
}

void ProcessBlockOnColumn( int size, block B1, block B2)
{
    for(int i=0; i < B2.size; i++) {
        for(int j=0; j < B1.height; j++) {
            A[B1.start+j*size+i] /= A[B2.start+i*size+i];
            for(int k = i+1; k < B2.size; k++) {
                A[B1.start+j*size+k] += -A[B1.start+j*size+i] * A[B2.start+i*size+k];
            }
        }
    }
}

void ProcessBlockOnRow( int size, block B1, block B2)
{
    for(int i=0; i < B2.size; i++)
        for(int j=i+1; j < B2.size; j++)
            for(int k=0; k < B1.size; k++)
                A[B1.start+j*size+k] += -A[B2.start+j*size+i] * A[B1.start+i*size+k];
}

void ProcessInnerBlock( int size, block B1, block B2, block B3)
{
    for(int i=0; i < B3.size; i++){
        for(int j=0; j < B1.height; j++){
            for(int k=0; k < B2.size; k++){
                A[B1.start+j*size+k] += -A[B3.start+j*size+i] * A[B2.start+i*size+k];
}

void checkResult( vector<double> &A2, int size) 
{
    int errors = 0;
    double temp2;
    double *L = new double[size*size];
    double *U = new double[size*size];
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
            if( (A2[i*size+j]-temp2) / A2[i*size+j] > 0.1 || (A2[i*size+j]-temp2) / A2[i*size+j] < -0.1 ){
                printf("error:[%d][%d]\n", i, j);
                printf("\t %f =/= %f \n", A2[i*size+j], temp2);
                errors++;
            }
        }
    if(errors > 0){
        printf("A:\n");
        Print_Matrix(A, size, size);
        printf("A2:\n");
        Print_Matrix(A2, size, size);
    }

    printf("Errors = %d \n", errors);
}

void Print_Matrix(vector<double> &v, int size1, int size)
{
    printf( "\n" );
    for(int i = 0; i < size1; i++){
        for(int j = 0; j < size; j++)
            printf( "%5.2f, ", v[i*size + j] );
        printf( "\n" );
    }
    printf( "\n" );
}

void InitMatrix3( int size)
{
    vector<future<void>> futures;
    futures.reserve(size);
    init_action innerLoop;
    for(int i = 0; i < size; i++)
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
    //for(long long i = 0; i < size; i++)
    //    for(long long j = 0; j < size; j++)
    //        for(long long k = 0; k < size; k++)
    //            A[i*size + j] += L[i*size + k] * U[k*size + j];
    for(int i = 0; i < size; i++) {
        futures.push_back( async(innerLoop, hpx::find_here(), i, size));
        //initLoop(i, size, L, U);
    }
    wait(futures);

}
void initLoop(int i, int size) {
    for(int j = 0; j < size; j++)
        for(int k = 0; k < size; k++)
            A[i*size + j] += L[i*size + k] * U[k*size + j];
}
