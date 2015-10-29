#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#define  MIN(x, y)  (x < y)?x:y
#define  MAX(x, y)  (x > y)?x:y

inline int similarity(char x, char y) {
    return((x==y)?2:-1);
}

unsigned long GetTickCount()
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec * 1000000) + (tv.tv_usec);
}

int find_array_max(int array[],int length);
char** openFile(char* fileName, unsigned int &N, unsigned int &num);
void sw( int **H, char *seq_a, char *seq_b, int chunk_size, int size_a, int size_b);

int ind, gap;

int main(int argc, char** argv)
{
    if(argc!=5){
        printf("Please provide correct number of input arguments: \n");
        printf("\t1:gap\n\t2:filename\n\tsequence_A\n\t3:filename sequence_B\n\t4:chunk_size\n");
        exit(0);
    }
    gap = atof(argv[1]);
    int chunk_size = atoi(argv[4]);
    unsigned int size_a, num_a, size_b, num_b;
    char **seq_a = openFile(argv[2], size_a, num_a);
    char **seq_b = openFile(argv[3], size_b, num_b);
    double t1,t2,total = 0;

    int **H = new int*[size_a + 1];
    for(int i=0; i< size_a + 1; i++)
        H[i] = new int[size_b + 1];

    for(int seq_num=0; seq_num < num_a; seq_num++) {
        for(int c=0; c<=size_a; c++){
            for(int d=0; d<=size_b; d++){
                H[c][d]=0;
            }
        }
        t1 = GetTickCount();
        sw(H, seq_a[seq_num], seq_b[seq_num], chunk_size, size_a, size_b);
        t2 = GetTickCount();
        total += t2 - t1; 
    }
    printf("Time for Smith-Watterman in secs: %f \n", total/1000000);
}

void sw(int **H, char *seq_a, char *seq_b, int chunk_size, int size_a, int size_b)
{
    int num_waves = size_a + size_b - 1; //note: was +1
    int elements, np, mp;
    int temp[4];

    for(int wave = 0; wave < num_waves; ++wave) {
        if(wave < size_a-1) {
            elements = wave+1;
            np = wave+1;
            mp = 1;
        } else if(wave < size_b) {
            elements = size_a;
            np = size_a;
            mp = wave-size_a+2;
        } else {
            elements = size_a+size_b-1-wave;
            np = size_a;
            mp = wave-size_a+2;
        }
        for(int i = 0; i < elements; i++) {
            temp[0] = H[(np-i)-1][(mp+i)-1] + similarity( seq_a[(np-i)-1], seq_b[(mp+i)-1] );
            temp[1] = H[(np-i)-1][(mp+i)  ] - gap;
            temp[2] = H[(np-i)  ][(mp+i)-1] - gap;
            temp[3] = 0;
            H[(np-i)][(mp+i)] = find_array_max( temp, 4 );
        }
    }
}

char** openFile(char* fileName, unsigned int &size, unsigned int &num)
{
    char **seq;
    FILE *fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Error: could not open %s!\n", fileName);
        exit(0);
    }
    fscanf(fp, "%d %d\n", &num, &size);
    seq = new char*[num];

    for(int a=0; a < num; a++) {
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

