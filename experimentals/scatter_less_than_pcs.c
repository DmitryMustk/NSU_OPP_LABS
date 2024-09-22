#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int* arr = malloc(sizeof(int) * 4);
	int* partArr = malloc(sizeof(int));

	if (rank == 0) {
		arr[0] = 1;
		arr[1] = 2;
		arr[3] = 3;
		arr[4] = 4;
	}
	MPI_Scatter(arr, 1, MPI_INT, partArr, 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("%d on %d process\n", *partArr, rank);

	MPI_Finalize();
}
