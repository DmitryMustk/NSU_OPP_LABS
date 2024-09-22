#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
	double* vec		   = malloc(sizeof(double) * 6);
	double* portionVec = malloc(sizeof(double) * 2);
	
	MPI_Init(NULL, NULL);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		for (int i = 0; i < 6; ++i) {
			vec[i] = 1.0;
		}
	}

	if (rank == 0) {
		MPI_Send(vec, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(portionVec, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("%f\n", portionVec[0]);
	}

	MPI_Finalize();
	return 0;
}
