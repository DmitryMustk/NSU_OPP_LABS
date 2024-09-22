#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void mulVecOnScalar(double* vec, int vecSize, double scalar, int rank) {
	#pragma omp parallel for 
	for (int i = 0; i < vecSize; ++i) {
		vec[i] *= scalar;
		printf("iter %d, from %d thread\n", i, rank);
	}
}

double getDotProduct(double* vec1, double* vec2, int vecSize) {
	double res = 0;
	#pragma omp parallel for reduction(+:res) 
	for (int i = 0; i < vecSize; ++i) {
		res += vec1[i] * vec2[i];
	}
	return res;
}

void fillVec(double* vec, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		vec[i] = 2.0;
	}
}

int main(void) {
	int vecSize = 12;
	double* vec1 = malloc(sizeof(double) * vecSize);
	double* vec2 = malloc(sizeof(double) * vecSize);
	
	fillVec(vec1, vecSize);
	fillVec(vec2, vecSize);

	#pragma omp parallel
	{
		int rank = omp_get_thread_num();
		double localDP = getDotProduct(vec1, vec2, vecSize);
		#pragma omp barrier
		#pragma omp single
		printf("res: %f\n", localDP);
	}
	return 0;
}

