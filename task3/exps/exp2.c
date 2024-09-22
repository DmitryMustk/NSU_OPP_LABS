#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void mulVecOnScalar(double* vec, int vecSize, double scalar, int rank) {
	#pragma omp for 
	for (int i = 0; i < vecSize; ++i) {
		vec[i] *= scalar;
		printf("iter %d, by %d thread\n", i, rank);
	} 
}

int main(void) {
	int vecSize = 12;
	double* vec = malloc(sizeof(double) * vecSize);
	double scalar = 2.0;
	#pragma omp parallel
	{
		int rank = omp_get_thread_num();
		mulVecOnScalar(vec, vecSize, scalar, rank);
	}
	return 0;
}
