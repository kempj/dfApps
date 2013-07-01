#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#define  MIN(x, y)  (x < y)?x:y
#define  MAX(x, y)  (x > y)?x:y

int similarity(char x, char y) {
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

int ind, gap;

int main(int argc, char** argv)
{
    if(argc!=5){
        printf("Please provide correct number of input arguments: \n");
        printf("1:gap 2:filename sequence_A  3:filename sequence_B 4:chunk_size\n");
        exit(0);
    }
    gap = atof(argv[1]);
    int chunk = atoi(argv[4]);
    unsigned int N_a, num_a, N_b, num_b;
    char **seq_a = openFile(argv[2], N_a, num_a);
    char **seq_b = openFile(argv[3], N_b, num_b);
    double t1,t2,total;

    int **H = new int*[N_a + 1];
    for(int i=0; i< N_a + 1; i++)
        H[i] = new int[N_b + 1];
    for(int i=0;i<=N_a;i++){
        for(int j=0;j<=N_b;j++){
            H[i][j]=0;
        }
    }

    int temp[4];
    int H_max = 0, i_max=0, j_max=0;
    int current_i, current_j, next_i, next_j;
    int tick=0;
    int waves = N_a + N_b +1; 
    int wave, elements, np, mp;
    int min = 0;

    for(int a=0; a < num_a; a++) {
        for(int c=0;c<=N_a;c++){
            for(int d=0;d<=N_b;d++){
                H[c][d]=0;
            }
        }
        t1 = omp_get_wtime();
#pragma omp parallel firstprivate(a, gap, waves) private(temp, wave, i) shared(np, mp, elements)
        {
#pragma omp master
            {
                for(wave = 0; wave < waves; ++wave) {
                    if(wave < N_a-1) {
                        elements = wave+1;
                        np = wave+1;
                        mp = 0+1;
                    } else if(wave < N_b) {
                        elements = N_a;
                        np = N_a-1+1;
                        mp = wave-(N_a-1)+1;
                    } else {
                        elements = N_a-1-(wave-N_b);
                        np = N_a-1+1;
                        mp = wave-(N_a-1)+1;
                    }
                    for(int ii = 0; ii < elements; ii+=chunk) {
                        min = MIN(elements,ii + chunk);
#pragma omp task firstprivate(ii, np, mp, chunk, elements)
                        {
                            for (int i = ii; i < min; i++) {
                                temp[0] = H[(np-i)-1][(mp+i)-1] + similarity(seq_a[a][(np-i)-1],seq_b[a][(mp+i)-1]); 
                                temp[1] = H[(np-i)-1][(mp+i)]-gap;                  
                                temp[2] = H[(np-i)][(mp+i)-1]-gap;
                                temp[3] = 0;
                                H[(np-i)][(mp+i)] = find_array_max(temp,4);
                            }
                        } 
                    } 
#pragma omp taskwait
                }
            }
        }
        t2 = omp_get_wtime();
        total += t2 - t1; 
    }
    printf("Time for Smith-Watterman in secs: %f \n", total);
}

char** openFile(char* fileName, unsigned int &N, unsigned int &num)
{
    char **seq;
    FILE *fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Error: could not open %s!\n", fileName);
        exit(0);
    }
    fscanf(fp, "%d %d\n", &num, &N);
    seq = new char*[num];

    for(int a=0; a < num; a++) {
        seq[a] = new char[N+1];
        fscanf(fp, "%s\n", seq[a]);
    }
    fclose(fp);
    return seq;
}

int find_array_max(int array[],int length)
{
    int max = array[0];
    ind = 0;
    for(int i = 1; i<length; i++){
        if(array[i] > max){
            max = array[i];
            ind = i; 
        }
    }
    return max;
}

