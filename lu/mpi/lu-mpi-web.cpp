#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <unistd.h>
#include <assert.h>
using namespace std;

class block {
    public:
        block(int H, int W);
        block() {}
        int height;
        int width;
        double *data;
        void print();
        void randInit();

        bool isReady();
        void wait();
        MPI_Request request;
        bool processed;
};

block::block(int H, int W) {
    height = H;
    width = W;
    data = new double[height * width];
}

void block::randInit() {
    for( int i = 0; i < height; i++) {
        for( int j = 0; j < width; j++) {
            data[i*width + j] = 5 * (rand()/(double)RAND_MAX);
        }
    }
}

void block::print() {
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            cout << data[i*width + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

bool block::isReady() {
    int flag = 0;
    MPI_Status status;
    MPI_Test(&request, &flag, &status);
    return (flag != 0);
}

void block::wait() {
    MPI_Wait(&request, NULL);
}

int getBlockDim(int rowOrCol, int size, int numBlocks) {
    int blockSize = size/numBlocks;
    int H;
    int numLargeBlocks = size % numBlocks;
    if( rowOrCol <= numLargeBlocks)
        H = blockSize + 1;
    else
        H = blockSize;
    return H;
}

block** init(int size, int numBlocks, int numLocalRows, int numNodes, int nodeID) {
    block** rowList = new block*[numLocalRows];
    int H, W;
    for(int i = 0; i < numLocalRows; i++) {
        int globalRow = i * numNodes + nodeID;
        H = getBlockDim(globalRow, size, numBlocks);
        rowList[i] = new block[numBlocks];
        for(int j = 0; j < numBlocks; j++) {
            W = getBlockDim(j, size, numBlocks);
            rowList[i][j] = block(H, W);
            rowList[i][j].randInit();
        }
    }
    return rowList;
}

void calcDiag(block &dBlock) {
    double *A = dBlock.data;
    int size = dBlock.width;
    for(int i = 0; i < size; i++) {
        for(int j = i+1; j < size; j++){
            A[j*size+i] /= A[i*size+i];
            for(int k = i+1; k < size; k++) {
                A[j*size+k] -= A[j*size+i] * A[i*size+k];
            }
        }
    }
}

void calcColBlock(block diag, block &previous) {
    double *A = previous.data;
    double *D = diag.data;
    for(int i=0; i < diag.width; i++) {
        for(int j=0; j < previous.height; j++) {
            A[j*previous.width+i] /= D[i*diag.width+i];
            for(int k = i+1; k < diag.width; k++) {
                A[j*previous.width+k] += -A[j*previous.width+i] * D[i*diag.width+k];
            }
        }
    }
}

void calcRowBlock(block diag, block &previous) {
    double *A = previous.data;
    double *D = diag.data;
    for(int i=0; i < diag.width; i++)
        for(int j=i+1; j < diag.width; j++)
            for(int k=0; k < previous.width; k++) {
                A[j*previous.width+k] += -D[j*diag.width+i] * A[i*previous.width+k];
            }
}

void calcInnerBlock(block row, block col, block &previous) {
    double *A = previous.data;
    double *R = row.data, *C = col.data;
    for(int i=0; i < row.width; i++)
        for(int j=0; j < previous.height; j++)
            for(int k=0; k < col.width; k++)
                A[j*previous.width+k] += -R[j*row.width+i] * C[i*col.width+k];
}

class receive {
    public:
        receive(int blocks, int nodes, int MatSize) :numBlocks(blocks), numNodes(nodes), size(MatSize){}
        block* operator()(int iter, int row, int col);
    private:
        int numBlocks, numNodes, size;
};

block* receive::operator()(int iter, int row, int col) {
    int height = getBlockDim(row, size, numBlocks);
    int width  = getBlockDim(col, size, numBlocks);
    block *B = new block(height, width);
    int size = B->height * B->width;
    MPI_Irecv(B->data, size, MPI_DOUBLE, row%numNodes, col + row*numBlocks + iter*numBlocks*numBlocks, MPI_COMM_WORLD, &B->request);
    B->processed = false;
    return B;
}

class broadcast {
    public:
        broadcast(int ID, int nNodes, int nBlocks) :nodeID(ID), numNodes(nNodes),  numBlocks(nBlocks){}
        void operator() (int iter, int row, int col, block B);
    private:
        int nodeID, numNodes, numBlocks;
};

void broadcast::operator()(int iter, int row, int col, block B){
    MPI_Request request;
    for(int i = 0; i < numNodes; i++) {
        if(i != nodeID){
            int size = B.height * B.width ;
            int tag = col + numBlocks * row + numBlocks * numBlocks * iter;
            cout << "node " << nodeID << " - block sent to node " << i << " of " << numNodes << endl;
            MPI_Isend(B.data, size, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &request);
        }
    }
}

void printList(block **localRowList, int nodeID, int localNumRows, int numBlocks) {
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(nodeID * 500000);
    for(int i = 0; i < localNumRows; i++) {
        for(int j = 0; j < numBlocks; j++) {
            cout << "node " << nodeID << endl << "Block " << i << ", " << j << endl;
            localRowList[i][j].print();
            cout << endl;
        }
    }
}

int main(int argc, char* argv[]){
    int numBlocks = 0, size = 0;
    int nodeID = 0, numNodes = 1;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);
    MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
    if(nodeID == 0) {
        if(argc > 2) {
            size = atoi(argv[1]);
            numBlocks = atoi(argv[2]);
            cout << "Size = " << size << ", number of blocks = " << numBlocks << endl;
        }
    }
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(size <= 0 || numBlocks <= 0) {
        cout << "input is:" << endl << "[size numBlocks]" << endl;
        MPI_Finalize();
            return 1;
    }
    if(numBlocks > 1290) {
        if(nodeID == 0) {
            cout << "The number of blocks is too large" << endl
                << "The current implementation only supports up to 1290" << endl;
        }
        MPI_Finalize();
        return 1;
    }
    if(numBlocks < numNodes) {
        if(nodeID == 0) 
            cout << "too many nodes for problem size" << endl;
        MPI_Finalize();
        return 1;
    }
    int localNumRows = numBlocks/numNodes;
    if(nodeID < numBlocks % numNodes) {
        localNumRows += 1;
    }
    block **localRowList = init(size, numBlocks, localNumRows, numNodes, nodeID);
    broadcast broadcastBlock(nodeID, numNodes, numBlocks);
    receive recvBlock(numBlocks, numNodes, size);

    block **topRow;
    block *diagBlock;

    for(int i = 0; i < numBlocks; i++){
        int localRowNum = i/numNodes;

        if(nodeID == i % numNodes) {
            cout << "node " << nodeID << " is calculating and sending row " << i << " of " << numBlocks <<  endl;
            
            calcDiag(localRowList[localRowNum][i]);

            if(numNodes > 1)
                broadcastBlock(i, i, i, localRowList[localRowNum][i]);

            for(int j=i+1; j < numBlocks; j++) {
                calcColBlock(localRowList[localRowNum][i], localRowList[localRowNum][j]);
                if(numNodes > 1)
                    broadcastBlock(i,i,j,localRowList[localRowNum][j]);
            }
            for(int j = localRowNum; j < localNumRows; j++) {
                calcRowBlock(localRowList[localRowNum][i], localRowList[j][i]);
            }
            for(int j = localRowNum+1; j < localNumRows; j++) {//tasks could split off here.
                for(int k = i+1; k < numBlocks; k++) {//or here
                    calcInnerBlock(localRowList[j][i], localRowList[localRowNum][k], localRowList[j][k]);
                }
            }
            cout << "node " << nodeID << " is the top row and is done calculating iteration " << i << endl;

        } else {
            cout << "node " << nodeID << " is receiving row " << i << endl;
            topRow = new block*[numBlocks];
            diagBlock = recvBlock(i,i,i);
            for(int j=i+1; j < numBlocks; j++) {
                topRow[j] = recvBlock(i,i,j);
            }

            cout << "node " << nodeID << " is waiting on diag " << i << endl;
            diagBlock->wait();
            cout << "node " << nodeID << " received diag: " << endl;

            for(int j = localRowNum; j < localNumRows; j++)
                calcRowBlock(*diagBlock, localRowList[j][i]);

            int numCompleted = 0;
            while( numCompleted < numBlocks - (i+1)) {
                for(int k = i+1; k < numBlocks; k++) {
                    if(topRow[k]->isReady() && topRow[k]->processed == false) {
                        for(int j = localRowNum; j < localNumRows; j++) {//tasks could split off here
                            calcInnerBlock(localRowList[j][i], *topRow[k], localRowList[j][k]);
                        } 
                        topRow[k]->processed = true;
                        numCompleted += 1;
                        cout << "node " << nodeID << " calculated column " << k << endl;
                    }
                }
            }
            delete [] topRow;
            cout << "node " << nodeID << " is done calculating iteration " << i << endl;
        }
    }
    MPI_Finalize();
    return 0;
}
