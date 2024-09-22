#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

uint64_t calculateS(double* vec1, double* vec2, size_t vecSize) {
	uint64_t s = 0;
	for (size_t i = 0; i < vecSize; ++i) {
		for (size_t j = 0; j < vecSize; ++j) {
			s += vec1[i] * vec2[j];
		}
	}
	return s;
}

void fillVecs(double* vec1, double* vec2, size_t vecSize) {
	for (size_t i = 0; i < vecSize; ++i) {
		vec1[i] = 1.0;
		vec2[i] = 2.0;
	}
}

int main(void) {
	size_t vecSize = 120000;
	double* vec1 = malloc(sizeof(double) * vecSize);
	double* vec2 = malloc(sizeof(double) * vecSize);

	fillVecs(vec1, vec2, vecSize);
	printf("Res: %ld\n", calculateS(vec1, vec2, vecSize));
	return 0;
}
