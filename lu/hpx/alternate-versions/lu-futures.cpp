// Static blocked LU Decomposition

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <future>

//#include "lu_utils.h"
#include "lu-local.h"
#include <vector>

using std::vector;
using std::future;
using std::async;

void LU( int size, int numBlocks);

void stage1( int size, int offset, vector<vector<block>> &blocks);
void stage2( int size, int offset, vector<vector<block>> &blocks);
void stage3( int size, int offset, vector<vector<block>> &blocks);

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
        getBlockList(blockList, numBlocks, size);

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

void stage1( int size, int offset, vector<vector<block>> &blocks)
{
    ProcessDiagonalBlock( size, blocks[offset][offset] );
}

void stage2( int size, int offset, vector<vector<block>> &blocks)
{
    int numBlocks = blocks[0].size();
    vector<future<block>> futures;

    for(int i=offset + 1; i<numBlocks; i++){
        futures.push_back( async( std::launch::async, ProcessBlockOnColumn, size, blocks[i][offset], blocks[offset][offset]));
        futures.push_back( async( std::launch::async, ProcessBlockOnRow, size, blocks[offset][i], blocks[offset][offset]));
    }
    for(int i = 0; i < futures.size(); i++)
        futures[i].wait();
}

void stage3( int size, int offset, vector<vector<block>> &blocks)
{
    int numBlocks = blocks[0].size();
    vector<future<block>> futures;
    futures.reserve(blocks.size() * blocks.size());
    for(int i=offset+1; i < numBlocks; i++) {
        for(int j=offset+1; j < numBlocks; j++){
            futures.push_back( async( std::launch::async, ProcessInnerBlock, size, blocks[i][j], blocks[offset][j], blocks[i][offset]));
        }
    }
    for(int i = 0; i < futures.size(); i++)
        futures[i].wait();
}
