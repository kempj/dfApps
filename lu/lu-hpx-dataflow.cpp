// Static blocked LU Decomposition

#include <stdio.h>
#include <hpx/hpx_init.hpp>
#include <hpx/include/threads.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrapped.hpp>

#include "lu-local.h"
#include <iostream>
using std::cout;
using std::endl;

using std::vector;
using hpx::util::unwrapped;
using hpx::lcos::shared_future;
using hpx::lcos::wait_all;
using hpx::async;
using hpx::lcos::local::dataflow;
using hpx::when_all;
using hpx::make_ready_future;
vector<double> A;

void LU( int size, int numBlocks)
{
    vector<vector<block>> blockList;
    vector<vector<vector<shared_future<block>>>> dfArray(2);
    shared_future<block> *diag_block, *first_col;
    shared_future<int> fsize = make_ready_future(size);

    auto rowOp = unwrapped(&ProcessBlockOnRow);
    auto colOp = unwrapped(&ProcessBlockOnColumn);
    auto innerOp = unwrapped(&ProcessInnerBlock);
    auto diagOp = unwrapped(&ProcessDiagonalBlock);

    getBlockList(blockList, numBlocks, size);

    for(int i = 0; i < 2; i++){
        dfArray[i].resize( numBlocks );
        for(int j = 0; j < numBlocks; j++){
            dfArray[i][j].resize( numBlocks, make_ready_future( block()));
        }
    }
    vector<vector<shared_future<block>>> &last_array = dfArray[0];
    vector<vector<shared_future<block>>> &current_array = dfArray[1];

    //converts blocks to shared_futures for dataflow input in dfArray[0]
    dfArray[0][0][0] = async( ProcessDiagonalBlock, size, blockList[0][0] );
    diag_block = &dfArray[0][0][0];
    for(int i = 1; i < numBlocks; i++) {
        dfArray[0][0][i] = dataflow( rowOp, fsize, make_ready_future( blockList[0][i] ), *diag_block);
    }
    for(int i = 1; i < numBlocks; i++) {
        dfArray[0][i][0] = dataflow( colOp, fsize, make_ready_future( blockList[i][0] ), *diag_block);
        first_col = &dfArray[0][i][0];
        for(int j = 1; j < numBlocks; j++) {
            dfArray[0][i][j] = dataflow( innerOp, fsize, make_ready_future( blockList[i][j]), dfArray[0][0][j], *first_col);
        }
    }
    //bulk of work, entirely in shared_futures
    for(int i = 1; i < numBlocks; i++) {
        last_array = dfArray[i-1];
        current_array = dfArray[i];

        current_array[i][i] = dataflow( diagOp, fsize, last_array[i][i]);
        diag_block = &current_array[i][i];
        for(int j = i + 1; j < numBlocks; j++){
            current_array[i][j] = dataflow( rowOp, fsize, last_array[i][j], *diag_block);
        }
        for(int j = i + 1; j < numBlocks; j++){
            current_array[j][i] = dataflow( colOp, fsize, last_array[j][i], *diag_block);
            first_col = &current_array[j][i];
            for(int k = i + 1; k < numBlocks; k++) {
                current_array[j][k] = dataflow( innerOp, fsize, last_array[j][k], current_array[i][k], *first_col);
           }
        }
    }
    wait_all(current_array[numBlocks-1][numBlocks-1]);
}

int hpx_main (int argc, char *argv[])
{
    vector<double> originalA;
    int size = 1000;
    int numBlocks = 10;
    unsigned long t1, t2;
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
    if(runCheck) {
        printf("Error checking enabled\n");
        originalA.reserve(size*size);
        for(int i = 0; i < size * size; i++) {
            originalA[i] = A[i];
        }
    }
    t1 = GetTickCount();
    if(numBlocks == 1) {
        ProcessDiagonalBlock( size, block(size, 0, size));
    } else if( numBlocks > 1) {
        LU( size, numBlocks);
    } else { 
        printf("Error: numBlocks must be greater than 0.\n");
    }
    t2 = GetTickCount();
    printf("Time for LU-decomposition in secs: %f \n", (t2-t1)/1000000.0);
    
    if(runCheck) {
        checkResult( originalA,  size );
    }
    return hpx::finalize();
}

int main(int argc, char *argv[])
{
    using namespace boost::assign;
    std::vector<std::string> cfg;
    cfg += "hpx.os_threads=" +
        boost::lexical_cast<std::string>(hpx::threads::hardware_concurrency());

//    printf("argc = %d\n", argc);
    return hpx::init(argc, argv, cfg);
}
