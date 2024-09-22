#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
	MPI_Init(NULL, NULL);

	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		printf("%d procs\n", size);
	}
	
	MPI_Finalize();
	return 0;
}
