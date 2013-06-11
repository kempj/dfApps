// Static blocked LU Decomposition

#include <stdio.h>
#include <stdlib.h>
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

struct block {
    int size;
    int start;
    block(int size, int startAddress) : size(size), start(startAddress){}
};
void LU(vector<double> &A, int size, int numBlocks);
void checkResult(vector<double> &A, vector<double> &A2, int size);

void stage1(vector<double> &A, int size, vector<block> blocks);
void stage2(vector<double> &A, int size, int offset, vector<vector<block>> blocks);
void stage3(vector<double> &A, int size, int offset, vector<vector<block>> blocks);

void ProcessDiagonalBlock(vector<double> &A, int size, block B);
void ProcessBlockOnColumn(vector<double> &A, int size, block B1, block B2);
void ProcessBlockOnRow(vector<double> &A, int size, block B1, block B2);
void ProcessInnerBlock(vector<double> &A, int size, block B1, block B2, block B3);

void getBlockList(vector<vector<block>> &blocks, int size, int numBlocks);

void Print_Matrix(vector<double> &v, int numBlocks, int size);
void InitMatrix2(vector<double> &A, int size);
void InitMatrix3(vector<double> &A, int size);

//HPX_PLAIN_ACTION(ProcessBlockOnColumn, column_action);
//HPX_PLAIN_ACTION(ProcessBlockOnRow, row_action);
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
   
    //Print_Matrix(A, size, size);
    checkResult( A, A2, size);
    return 0;
}

void LU(vector<double> &A, int size, int numBlocks)
{
    double t1, t2;
    int blockSize , start;
    vector<vector<block>> blockList;

    for(int i=0; i < numBlocks; i++) 
        blockList.push_back(vector<block>());

    t1 = GetTickCount();
    if(numBlocks > 1) {
        getBlockList(blockList, size, numBlocks);

        for(int i = 0; i < blockList.size(); i++) {
            stage1( A, size, blockList.back() );
            stage2( A, size, i, blockList );
            stage3( A, size, i, blockList );
        }
        ProcessDiagonalBlock(A, size, blockList.back().back());
    }
    ProcessDiagonalBlock(A, size, block(size, 0));
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000);
}

void getBlockList(vector<vector<block>> &blockList, int size, int numBlocks)
{
    int blockSize, start;

    for(int i=0; i < numBlocks; i++) {
        if(i < size % numBlocks) {
            blockSize = size/numBlocks+1;
            start = (size/numBlocks+1)*i; //blockList.back().push_back( block( remain/numBlocks+1, ((remain/numBlocks+1)*i) + offset * size + offset));
        } else {
            blockSize = size/numBlocks;
            start = (size/numBlocks+1)*(size%numBlocks) + (size/numBlocks)*(i-size/numBlocks);
        }
        for(int j = 0; j < numBlocks; j++) {
            blockList[j].push_back( block( blockSize, start));
        }
                //remain/numBlocks, ((remain/numBlocks+1)*(remain%numBlocks) + (remain/numBlocks)*(i-remain%numBlocks) ) + offset * size + offset));
    }
}

void stage1(vector<double> &A, int size, vector<block> blocks)
{
    ProcessDiagonalBlock( A, size, blocks[0] );
}

void ProcessDiagonalBlock(vector<double> &A, int size, block B)
{
    for(int i = 0; i < B.size; i++)
        for(int j = i+1; j < B.size; j++){
            A[B.start+j*size+i] /= A[B.start+i*size+i];
            for(int k = i+1; k < B.size; k++)
                A[B.start+j*size+k] -= A[B.start+j*size+i] * A[B.start+i*size+k];
        }
}

void stage2(vector<double> &A, int size, int offset, vector<vector<block>> blocks)
{
    //column_action blockCol;
    //row_action blockRow;
    //vector<future<void>> futures;
    //futures.reserve(numBlocks*2);
//    int numBlocks = blocks[offset].size() - offset;
    int numBlocks = blocks[0].size();
    for(int i=1; i<numBlocks; i++){
        //futures.push_back( async<blockCol>( hpx::find_here(), A, 
        //                          (offset+start[i])*size + offset, 
        //                          offset*size + offset,
        //                          L1, L2, size ));
        ProcessBlockOnColumn( A, size, blocks[offset + i][offset], blocks[offset][offset]);
        ProcessBlockOnRow(    A, size, blocks[offset][offset + i], blocks[offset][offset]);
    }
    //wait(futures);
}

void ProcessBlockOnColumn(vector<double> &A, int size, block B1, block B2)
{
    for(int i=0; i < B1.size; i++)
        for(int j=0; j < B2.size; j++){
            A[B1.start+j*size+i] /= A[B2.start+i*size+i];
            for(int k = i+1; k < B1.size; k++)
                A[B1.start+j*size+k] += -A[B1.start+j*size+i] * A[B2.start+i*size+k];
        }
}

void ProcessBlockOnRow(vector<double> &A, int size, block B1, block B2)
{
    for(int i=0; i < B1.size; i++)
        for(int j=i+1; j < B1.size; j++)
            for(int k=0; k < B2.size; k++)
                A[B1.start + j*size+k] += -A[B2.start + j*size+i] * A[B1.start+i*size+k];
}

void stage3(vector<double> &A, int size, int offset, vector<vector<block>> blocks)
{
    int numBlocks = blocks[0].size();//offset].size() - offset;
    for(int i=offset+1; i < numBlocks; i++)
        for(int j=offset+1; j < numBlocks; j++){
            ProcessInnerBlock( A, size, 
                               blocks[i][j], 
                               blocks[offset][j], 
                               blocks[i][offset]);
        }
}

void ProcessInnerBlock(vector<double> &A, int size, block B1, block B2, block B3)
{
    for(int i=0; i < B1.size; i++)
        for(int j=0; j < B2.size; j++)
            for(int k=0; k < B3.size; k++)
                A[B1.start+j*size+k] += -A[B3.start+j*size+i] * A[B2.start+i*size+k];
}

void checkResult(vector<double> &A, vector<double> &A2, int size) 
{
    int temp2, errors = 0;
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
            if( (A2[i*size+j]-temp2) / A2[i*size+j] > 0.1 || (A2[i*size+j]-temp2) / A2[i*size+j] < -0.1 )
                errors++;
        }
    printf("Errors = %d \n", errors);
}

void Print_Matrix(vector<double> &v, int numBlocks, int size)
{
    printf( "\n" );
    for(int i = 0; i < numBlocks; i++){
        for(int j = 0; j < size; j++)
            printf( "%4.0f,", v[i*size + j] );
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
