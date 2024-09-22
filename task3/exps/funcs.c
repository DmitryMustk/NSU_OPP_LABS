#include <omp.h>
#include <stdio.h>

int main(void) {
	#pragma omp parallel 
	{
		int size = omp_get_num_threads();
		int rank = omp_get_thread_num();

		printf("hello from %d thread\n", rank);
	}
	return 0;

}
