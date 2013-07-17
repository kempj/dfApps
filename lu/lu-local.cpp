#include "lu_utils.h"

extern std::vector<double> A;

block ProcessDiagonalBlock( int size,  block B)
{
    for(int i = 0; i < B.size; i++) {
        for(int j = i+1; j < B.size; j++){
            A[B.start+j*size+i] /= A[B.start+i*size+i];
            for(int k = i+1; k < B.size; k++) {
                A[B.start+j*size+k] -= A[B.start+j*size+i] * A[B.start+i*size+k];
            }
        }
    }
    return  B;
}

block ProcessBlockOnColumn( int size, block B1, block B2)
{
    for(int i=0; i < B2.size; i++) {
        for(int j=0; j < B1.height; j++) {
            A[B1.start+j*size+i] /= A[B2.start+i*size+i];
            for(int k = i+1; k < B2.size; k++) {
                A[B1.start+j*size+k] += -A[B1.start+j*size+i] * A[B2.start+i*size+k];
            }
        }
    }
    return B1;
}

block ProcessBlockOnRow( int size, block B1, block B2)
{
    for(int i=0; i < B2.size; i++)
        for(int j=i+1; j < B2.size; j++)
            for(int k=0; k < B1.size; k++)
                A[B1.start+j*size+k] += -A[B2.start+j*size+i] * A[B1.start+i*size+k];
    return B1;
}

block ProcessInnerBlock( int size, block B1, block B2, block B3)
{
    for(int i=0; i < B3.size; i++)
        for(int j=0; j < B1.height; j++)
            for(int k=0; k < B2.size; k++)
                A[B1.start+j*size+k] += -A[B3.start+j*size+i] * A[B2.start+i*size+k];
    return B1;
}

