#include "lu_utils.h"

#include <stdio.h>

#include <hpx/hpx_init.hpp>
#include <hpx/include/threads.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/actions.hpp>

#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrapped.hpp>
#include <vector>

using hpx::util::unwrapped;
using std::vector;
using hpx::lcos::future;
using hpx::lcos::wait;
using hpx::async;
using hpx::lcos::local::dataflow;
using hpx::when_all;
using hpx::make_ready_future;

void LU( int size, int numBlocks);
void checkResult( vector<double> &originalA);
//block ProcessDiagonalBlock( int size, block B);
//block ProcessBlockOnColumn( int size, block B1, block B2);
//block ProcessBlockOnRow( int size, block B1, block B2);
//block ProcessInnerBlock( int size, block B1, block B2, block B3);
void getBlockList(vector<vector<block>> &blocks, int numBlocks, int size);
