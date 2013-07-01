#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char *argv[]) 
{
	if(argc < 5) {
		printf("Error: Too few arguments!\n");
		exit(0);
	}

	unsigned int i, j;
	FILE *fp;
    char* filename = argv[1];
    unsigned int num_alphabets = atoi(argv[2]);
    unsigned int num_sequences = atoi(argv[3]);
    unsigned int seq_size = atoi(argv[4]);
	srand(time(NULL));

	if(!(fp = fopen(argv[1], "w"))) {
		printf("Error: Could not open file!\n");
		exit(0);
	}
	
	fprintf(fp, "%d %d\n", num_sequences, seq_size);
	for(i = 0; i < num_sequences; ++i) {
		for(j = 0; j < seq_size; ++j) {
			fprintf(fp, "%c", (char) 65 + rand() % num_alphabets);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	return 0;
}
