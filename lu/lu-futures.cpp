// Static blocked LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <future>
#include <vector>

#include "lu_utils.h"

using std::vector;
using std::future;
using std::async;

void LU( int size, int numBlocks);
void checkResult( vector<double> &A2, int size);

void stage1( int size, int offset, vector<vector<block>> &blocks);
void stage2( int size, int offset, vector<vector<block>> &blocks);
void stage3( int size, int offset, vector<vector<block>> &blocks);

void ProcessDiagonalBlock( int size, block B);
void ProcessBlockOnColumn( int size, block B1, block B2);
void ProcessBlockOnRow( int size, block B1, block B2);
void ProcessInnerBlock( int size, block B1, block B2, block B3);

void getBlockList(vector<vector<block>> &blocks, int size, int numBlocks);

void Print_Matrix(vector<double> &v, int numBlocks, int size);
void InitMatrix3( int size);
void initLoop(int i, int size);

vector<double> A;

int main (int argc, char *argv[])
{
    unsigned long t1, t2;
    int size = 100, numBlocks = 1;
    bool runCheck = false;
    if( argc > 1 )
        size = atoi(argv[1]);
    if( argc > 2 )
        numBlocks = atoi(argv[2]);
    if( argc > 3 )
        runCheck = true;
    printf("size = %d, numBlocks = %d\n", size, numBlocks);
    A.resize(size*size, 0);
    //L.resize(size*size, 0);
    //U.resize(size*size, 0);
    InitMatrix3( size );
    vector<double> A2;
    A2.reserve(size*size);
    for(int i = 0; i < size * size; i++)
        A2[i] = A[i];
    LU( size, numBlocks);
   
    if(runCheck)
        checkResult( A2, size);
    return 0;
}

void LU( int size, int numBlocks)
{
    unsigned long t1, t2;
    vector<vector<block>> blockList;

    t1 = GetTickCount();
    if(numBlocks > 1) {
        getBlockList(blockList, size, numBlocks);

        for(int i = 0; i < blockList.size(); i++) {
            stage1( size, i, blockList );
            stage2( size, i, blockList );
            stage3( size, i, blockList );
        }
    } else {
        ProcessDiagonalBlock(size, block(size, 0, size));
    }
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000.0);
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

void stage1( int size, int offset, vector<vector<block>> &blocks)
{
    ProcessDiagonalBlock( size, blocks[offset][offset] );
}

void ProcessDiagonalBlock( int size, block B)
{
    for(int i = 0; i < B.size; i++)
        for(int j = i+1; j < B.size; j++){
            A[B.start+j*size+i] /= A[B.start+i*size+i];
            for(int k = i+1; k < B.size; k++)
                A[B.start+j*size+k] -= A[B.start+j*size+i] * A[B.start+i*size+k];
        }
}

void stage2( int size, int offset, vector<vector<block>> &blocks)
{
    int numBlocks = blocks[0].size();
    vector<future<void>> futures;

    for(int i=offset + 1; i<numBlocks; i++){
        futures.push_back( async(std::launch::async,  ProcessBlockOnColumn, size, blocks[i][offset], blocks[offset][offset]));
        futures.push_back( async(std::launch::async, ProcessBlockOnRow, size, blocks[offset][i], blocks[offset][offset]));
    }
    for(int i = 0; i < futures.size(); i++)
        futures[i].wait();
}

void ProcessBlockOnColumn( int size, block B1, block B2)
{
    for(int i=0; i < B2.size; i++) {
        for(int j=0; j < B1.height; j++){
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

void stage3( int size, int offset, vector<vector<block>> &blocks)
{
    int numBlocks = blocks[0].size();
    vector<future<void>> futures;
    futures.reserve(blocks.size() * blocks.size());
    for(int i=offset+1; i < numBlocks; i++) {
        for(int j=offset+1; j < numBlocks; j++){
            futures.push_back( async(std::launch::async, ProcessInnerBlock, size, blocks[i][j], blocks[offset][j], blocks[i][offset]));
        }
    }
    for(int i = 0; i < futures.size(); i++)
        futures[i].wait();
}

void ProcessInnerBlock( int size, block B1, block B2, block B3)
{
    for(int i=0; i < B3.size; i++){
        for(int j=0; j < B1.height; j++){
            for(int k=0; k < B2.size; k++){
                A[B1.start+j*size+k] += -A[B3.start+j*size+i] * A[B2.start+i*size+k];
            }
        }
    }
}
