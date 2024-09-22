#include <mpi/mpi.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

void fillVecs(double* vec1, double* vec2, size_t vecSize) {
	for (size_t i = 0; i < vecSize; ++i) {
		vec1[i] = 1.0;
		vec2[i] = 2.0;
	}
}

void fillVec(double* vec, size_t vecSize) {
	for (size_t i = 0; i < vecSize; ++i) {
		vec[i] = 1.0;
	}
}

void printVec(double* vec, size_t vecSize, int rank) {
	printf("(");
	for (size_t i = 0; i < vecSize; ++i) {
		printf("%f, ", vec[i]);
	}
	printf(") -- from %d process\n", rank);
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
	double* portionVec2 = malloc(sizeof(double) * portionSize);

	uint64_t s;

	if (rank == 0) {
		double* vec2 = malloc(sizeof(double) * vecSize);	
		fillVecs(vec1, vec2, vecSize);	
		for (int i = 1; i < size; ++i) {
			MPI_Send(vec1, vecSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(vec2 + i * portionSize, portionSize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
		}
		uint64_t res = calculateS(vec1, vecSize, vec2, portionSize);
		for (int i = 1; i < size; ++i) {
			uint64_t fromS = 0;
			MPI_Recv(&fromS, 1, MPI_UINT64_T, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			res += fromS;
		}
		printf("Res = %ld\n", res);
	}
	else {
		MPI_Recv(vec1, vecSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		MPI_Recv(portionVec2, portionSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		uint64_t s = calculateS(vec1, vecSize, portionVec2, portionSize);
		MPI_Send(&s, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	}
	MPI_Finalize();

	return 0;
}
