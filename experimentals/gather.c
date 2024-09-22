#include <mpi/mpi.h>
#include <stdlib.h>
#include <stdio.h>

void fillVec(int* vec, int vecSize) {
	for (int i = 0; i < 2; ++i) {
		vec[i] = 1;
	}
}

void printVec(int* vec, int vecSize, int rank) {
	printf("(");
	for (int i = 0; i < vecSize; ++i) {
		printf("%d, ", vec[i]);
	}
	printf(") -- from %d process\n", rank);
}


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int vecSize = 10;

	int* tmpVec = malloc(sizeof(int) * 2);
	int* tmp2Vec = malloc(sizeof(int) * vecSize);
	fillVec(tmpVec, vecSize);

	//MPI_Allgather(tmpVec, 2, MPI_INT, tmp2Vec, 2, MPI_INT, MPI_COMM_WORLD);
	//MPI_Gather(tmpVec, 2, MPI_INT, tmp2Vec, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Allgather(tmpVec, 2, MPI_INT, tmp2Vec, 2, MPI_INT, MPI_COMM_WORLD);
	
	printVec(tmp2Vec, vecSize, rank);
	MPI_Finalize();

	free(tmpVec);
	free(tmp2Vec);
	return 0;
}
