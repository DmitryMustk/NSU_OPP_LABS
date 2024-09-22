#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void fillVec(double* vec, size_t vecSize) {
	for (size_t i = 0; i < vecSize; ++i) {
		vec[i] = i + 1;
	}
}

void fillVecs(double* vec1, double* vec2, size_t vecSize) {
	for (size_t i = 0; i < vecSize; ++i) {
		vec1[i] = 1.0;
		vec2[i] = 2.0;
	}
}

int main(void) {
	size_t vecSize = 10;

	double* vec1 = malloc(sizeof(double) * vecSize);
	double* vec2 = malloc(sizeof(double) * vecSize);


	fillVecs(vec1, vec2, vecSize);
	//uint64_t s = computeS(vec1, vec2, vecSize);
	//uint64_t local_s = computeS(vec1, vec2, vecSize);
	//
	double res = 0;
	#pragma omp parallel
	{
		double local_s = 0;
		#pragma omp for reduce(+:res)
		for (size_t i = 0; i < vecSize; ++i) {
			local_s += vec1[i] * vec2[i];
		}
	}
	printf("res: %f\n", res);


	free(vec1);
	free(vec2);
	return 0;
}
