#include <mpi/mpi.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

void fillVecs(double* vec1, double* vec2, size_t vecSize) {
	for (size_t i = 0; i < vecSize; ++i) {
		vec1[i] = 1.0;
		vec2[i] = 2.0;
	}
}

uint64_t calculateS(double* vec1, size_t vec1Size, double* vec2, size_t vec2Size) {
	uint64_t s = 0;
	for (size_t i = 0; i < vec1Size; ++i) {
		for (size_t j = 0; j < vec2Size; ++j) {
			s += vec1[i] * vec2[j];
		}
	}
	return s;
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	size_t vecSize = 120000;
	size_t portionSize = vecSize / size;
	double* vec1 = malloc(sizeof(double) * vecSize);
	double* vec2 = malloc(sizeof(double) * vecSize);
	double* portionVec2 = malloc(sizeof(double) * portionSize);
	if (rank == 0) {
		fillVecs(vec1, vec2, vecSize);
	}
	MPI_Bcast(vec1, vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(vec2, portionSize, MPI_DOUBLE, portionVec2, portionSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	uint64_t res = 0;
	uint64_t s = calculateS(vec1, vecSize, portionVec2, portionSize);
	MPI_Reduce(&s, &res, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("Res: %ld\n", res);

	}

	MPI_Finalize();
	return 0;
}
