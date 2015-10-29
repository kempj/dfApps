#include <sw_hpx.h>

#include <hpx/hpx_init.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/actions.hpp>

#include <boost/program_options.hpp>
#include <boost/assign.hpp>

//#define  MIN(x, y)  (x < y)?x:y
//#define  MAX(x, y)  (x > y)?x:y

using std::vector;
using hpx::lcos::future;
using hpx::lcos::wait_all;
using hpx::async;

inline int similarity(char x, char y) {
    return((x==y)?2:-1);
}

unsigned long GetTickCount()
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec * 1000000) + (tv.tv_usec);
}

int gap, chunk_size;
int **H;

int hpx_main(int argc, char **argv)
{
    gap = atof(argv[1]);
    chunk_size = atoi(argv[4]);
    size_t size_a, num_a, size_b, num_b;
    char **seq_a = openFile(argv[2], size_a, num_a);
    char **seq_b = openFile(argv[3], size_b, num_b);
    double t1, t2, total = 0;

    H = new int*[size_a + 1];
    for(size_t i=0; i< size_a + 1; i++) {
        H[i] = new int[size_b + 1];
    }
    for(size_t seq_num=0; seq_num < num_a; seq_num++) {
        for(size_t i=0; i<=size_a; i++) {
            std::fill(H[i], H[i] + (size_b+1), 0);
        }
        t1 = GetTickCount();
        sw(seq_a[seq_num], seq_b[seq_num], size_a, size_b);
        t2 = GetTickCount();
        total += t2 - t1; 
    }
    printf("Time for Smith-Watterman in secs: %f \n", total / 1000000);

    return hpx::finalize();
}


int main(int argc, char** argv)
{
    if(argc!=5){
        printf("Please provide correct number of ia_offsetut arguments: \n");
        printf("\t1:gap\n\t2:filename\n\tsequence_A\n\t3:filename sequence_B\n\t4:chunk_size\n");
        exit(0);
    }

    using namespace boost::assign;
    std::vector<std::string> cfg;
    cfg += "hpx.os_threads=" +
        boost::lexical_cast<std::string>(hpx::threads::hardware_concurrency());

    return hpx::init(argc, argv, cfg);
}

void sw(char *seq_a, char *seq_b, int size_a, int size_b)
{
    int num_waves = size_a + size_b - 1;
    int elements, a_offset, b_offset, end;
    vector< future<void> > futures;

    for(int wave = 0; wave < num_waves; ++wave) {
        if(wave < size_a-1) {
            elements = wave+1;
            a_offset = wave+1;
            b_offset = 1;
        } else if(wave < size_b) {
            elements = size_a;
            a_offset = size_a;
            b_offset = wave-size_a+2;
        } else {
            elements = size_a+size_b-1-wave;
            a_offset = size_a;
            b_offset = wave-size_a+2;
        }
        for(int ii = 0; ii < elements; ii+=chunk_size) {
            end = MIN(elements,ii + chunk_size);
            futures.push_back(async(innerloop, ii, end, seq_a, seq_b, b_offset, a_offset));
        }
        wait_all(futures);
    }
}

void innerloop(int start, int end, char *seq_a, char *seq_b, int b_offset, int a_offset)
{
    int temp[4];
    for (int i = start; i < end; i++) {
        temp[0] = H[(a_offset-i)-1][(b_offset+i)-1] + similarity( seq_a[(a_offset-i)-1], seq_b[(b_offset+i)-1] );
        temp[1] = H[(a_offset-i)-1][(b_offset+i)  ] - gap;
        temp[2] = H[(a_offset-i)  ][(b_offset+i)-1] - gap;
        temp[3] = 0;
        H[a_offset-i][b_offset+i] = find_array_max( temp, 4 );
    }
} 

char** openFile(char* fileName, size_t &size, size_t &num)
{
    char **seq;
    FILE *fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Error: could not open %s!\n", fileName);
        exit(0);
    }
    fscanf(fp, "%zu %zu\n", &num, &size);
    seq = new char*[num];

    for(size_t a=0; a < num; a++) {
        seq[a] = new char[size+1];
        fscanf(fp, "%s\n", seq[a]);
    }
    fclose(fp);
    return seq;
}

inline int find_array_max(int array[],int length)
{
    int max = array[0];
    for(int i = 1; i<length; i++){
        if(array[i] > max){
            max = array[i];
        }
    }
    return max;
}

