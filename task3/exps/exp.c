#include <stdio.h>
#include <omp.h>

#define N 10000

float x[N];

int main(void) {
	printf("Point 0\n");
	#pragma omp parallel
	{
		int rank = omp_get_thread_num();
		int size = omp_get_num_threads();
		printf("hello from %d of %d\n", rank, size);
	}
	printf("Point 1\n");
	#pragma omp parallel 
	{
		printf("Hi!\n");
	}
	printf("Point 2\n");
	return 0;
}
