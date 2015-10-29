#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#define  MIN(x, y)  (x < y)?x:y
#define  MAX(x, y)  (x > y)?x:y



int ind;
int indx;
int gap;

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
void check_dep_condition(int task_id, int task_cnt_prev1, int task_cnt_prev2, int elements, int line, int chunk);

int main(int argc, char** argv)
{
    if(argc!=5){
        printf("Please provide correct number of input arguments: \n");
        printf("1:gap 2:filename sequence_A  3:filename sequence_B 4:chunk_size\n");
        exit(0);
    }
    gap = atof(argv[1]); // loading gap
    int i,ii,j,a,b,c,d;
    double t1,t2,total;
    int chunk; 
    chunk = atoi(argv[4]);

    //load sequences
    FILE *fpa;
    fpa = fopen(argv[2], "r");
    if (fpa == NULL) {
        printf("Error: could not open seqA file!\n");
        exit(0);
    }
    unsigned int N_a, num_a;
    fscanf(fpa, "%d %d\n", &num_a, &N_a);
    char **seq_a = malloc(num_a * sizeof(char*));
    for(a=0; a < num_a; a++) {
        seq_a[a] = malloc((N_a+1) * sizeof(char));
        fscanf(fpa, "%s\n", seq_a[a]);
    }
    fclose(fpa);

    FILE *fpb;
    fpb = fopen(argv[3], "r");
    if (fpb == NULL) {
        printf("Error: could not open seqB file!\n");
        exit(0);
    }

    unsigned int N_b, num_b;
    fscanf(fpb, "%d %d\n", &num_b, &N_b);
    char **seq_b = malloc(num_b * sizeof(char*));
    for(a=0; a < num_b; a++) {
        seq_b[a] = malloc((N_b+1) * sizeof(char));
        fscanf(fpb, "%s\n", seq_b[a]);
    }
    fclose(fpb);

    //initialize H
    int **H = malloc((N_a + 1) * sizeof(int *));

    for(i=0; i< N_a + 1; i++) {
        H[i] = malloc((N_b + 1) * sizeof(int));
    }
    for(i=0;i<=N_a;i++) {
        for(j=0;j<=N_b;j++) {
            H[i][j]=0;
        }
    }



    int temp[4];
    int wave_cell = ceil(N_a/chunk);
    int H_max = 0.;
    int i_max=0,j_max=0;
    int current_i,current_j;
    int next_i;
    int next_j;
    int tick=0;
    int line;
    int task_count = 0;
    int task_cnt_prev1, task_cnt_prev2;

    int waves = N_a + N_b +1; 
    int wave, elements, np, mp;
    int min = 0;
    int task_id = 0;
    int **tab;
    tab = malloc(waves * sizeof(int *));
    for(i=0; i<waves; i++) {
        tab[i] = malloc(wave_cell*sizeof(int));
    }
    for(i=0;i< waves;i++) {
        for(j=0;j< wave_cell ;j++) {
            tab[i][j]=0;
        }
    }
    for(a=0; a < num_a; a++) {
        for(c=0;c<=N_a;c++) {
            for(d=0;d<=N_b;d++) {
                H[c][d]=0;
            }
        }
        t1 = omp_get_wtime();
#pragma omp parallel firstprivate(a, gap, waves) private(temp, wave, ii, i, task_count) shared(np, mp, elements, line)
        {
#pragma omp master
            {
                for(wave = 0; wave < waves; ++wave) {
                    if(wave < N_a-1) {
                        elements = wave+1;
                        np = wave+1;
                        mp = 0+1;
                        line = -1;
                    } else if(wave < N_b) {
                        elements = N_a;
                        np = N_a-1+1;
                        mp = wave-(N_a-1)+1;
                        line = 0;
                    } else {
                        elements = N_a-1-(wave-N_b);
                        np = N_a-1+1;
                        mp = wave-(N_a-1)+1;
                        line = 1;
                    }

                    task_cnt_prev2 = task_cnt_prev1;
                    task_cnt_prev1 = task_count; 
                    task_count=0;

                    for(ii = 0; ii < elements; ii+=chunk) {
                        min = MIN(elements,ii + chunk);
                        task_id++;
                        if (line == -1 && elements <= chunk) {






#pragma omp task firstprivate(ii, np, mp, chunk, elements, min, task_id) \
                            in(task_id-1) in(task_id -2) \
                            out(task_id)
                            {
                                for (i = ii; i < min; i++) {
                                    temp[0] = H[(np-i)-1][(mp+i)-1] + similarity(seq_a[a][(np-i)-1],seq_b[a][(mp+i)-1]); 
                                    temp[1] = H[(np-i)-1][(mp+i)]-gap;                  
                                    temp[2] = H[(np-i)][(mp+i)-1]-gap;
                                    temp[3] = 0;
                                    H[(np-i)][(mp+i)] = find_array_max(temp,4);
                                }
                            } // task
                        } else if( line <= 0 && elements > chunk) {
#pragma omp task firstprivate(ii, np, mp, chunk, elements, min, task_id, task_cnt_prev1, task_cnt_prev2) \
                            in(task_id - task_cnt_prev1    ) \
                            in(task_id - task_cnt_prev1 - 1) \
                            in(task_id -(task_cnt_prev1 + task_cnt_prev2    )) \
                            in(task_id -(task_cnt_prev1 + task_cnt_prev2 - 1)) \
                            out(task_id)
                            {
                                for (i = ii; i < min; i++) {
                                    temp[0] = H[(np-i)-1][(mp+i)-1] + similarity(seq_a[a][(np-i)-1],seq_b[a][(mp+i)-1]); 
                                    temp[1] = H[(np-i)-1][(mp+i)]-gap;                  
                                    temp[2] = H[(np-i)][(mp+i)-1]-gap;
                                    temp[3] = 0;
                                    H[(np-i)][(mp+i)] = find_array_max(temp,4);
                                }
                            } // task
                        } else if( line > 0 && elements >= chunk) {
#pragma omp task firstprivate(ii, np, mp, chunk, elements, min, task_id, task_cnt_prev1, task_cnt_prev2) \
                            in(task_id-task_cnt_prev1) \
                            in(task_id-task_cnt_prev1+1) \
                            in(task_id-(task_cnt_prev1 + task_cnt_prev2)) \
                            in(task_id-(task_cnt_prev1 + task_cnt_prev2 +1)) \
                            out(task_id)
                            {
                                for (i = ii; i < min; i++) {
                                    temp[0] = H[(np-i)-1][(mp+i)-1] + similarity(seq_a[a][(np-i)-1],seq_b[a][(mp+i)-1]); 
                                    temp[1] = H[(np-i)-1][(mp+i)]-gap;                  
                                    temp[2] = H[(np-i)][(mp+i)-1]-gap;
                                    temp[3] = 0;
                                    H[(np-i)][(mp+i)] = find_array_max(temp,4);
                                }
                            }
                        } else {
#pragma omp task firstprivate(ii, np, mp, chunk, elements, min, task_id) \
                            in(task_id-1) in(task_id-2) in(task_id-3) \
                            out(task_id)
                            {
                                for (i = ii; i < min; i++) {
                                    temp[0] = H[(np-i)-1][(mp+i)-1] + similarity(seq_a[a][(np-i)-1],seq_b[a][(mp+i)-1]); 
                                    temp[1] = H[(np-i)-1][(mp+i)]-gap;                  
                                    temp[2] = H[(np-i)][(mp+i)-1]-gap;
                                    temp[3] = 0;
                                    H[(np-i)][(mp+i)] = find_array_max(temp,4);
                                }
                            }
                        }
                        task_count++;
                    } //for loop
                }
            }
        }
        t2 = omp_get_wtime();
        total += t2 - t1; 
    }
    printf("Time for Smith-Watterman in secs: %f \n", total);
} // END of main

void check_dep_condition(int task_id, int task_cnt_prev1, int task_cnt_prev2, int elements, int line, int chunk)
{
    printf("line=%d, elements=%d\n",line,elements);
    if (line == -1 && elements <= chunk)
        indx = 1;
    if (line <= 0 && elements > chunk)
        indx = 2;
    if (line > 0 && elements >= chunk)
        indx = 3;
    if (line > 0 && elements < chunk)
        indx = 4;
    else
        indx = 0;
}

int find_array_max(int array[],int length)
{
    int max = array[0];            // start with max = first element
    ind = 0;
    int i;

    for(i = 1; i<length; i++){
        if(array[i] > max){
            max = array[i];
            ind = i; 
        }
    }
    // printf("max=%d\n", max);
    // printf(".........\n");
    return max;                    // return highest value in array
}

