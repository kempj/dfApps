#include <vector>

extern std::vector<double> A;

struct block {
    int size;
    int start;
    int height;
    block(int size, int startAddress, int H) : size(size), start(startAddress), height(H){}
    block() : size(0), start(0), height(0){}
};

unsigned long GetTickCount();
void InitMatrix3( std::vector<double> &L, std::vector<double> &U, int size);
void initLoop( std::vector<double> const& L, std::vector<double> const& U, int i, int size);
void Print_Matrix( std::vector<double> &v, int size);
void checkResult( std::vector<double> &originalA, int size );

