#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_OF_DIMS 2

int main(void) {
	MPI_Init(NULL, NULL);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int dims[] = {2, 3};
	int periods[] = {0, 0};
	int cords[NUM_OF_DIMS];

	MPI_Comm MPI_COMM_GRID;
	MPI_Cart_create(MPI_COMM_WORLD, NUM_OF_DIMS, dims, periods, 0, &MPI_COMM_GRID);
	MPI_Cart_coords(MPI_COMM_GRID, rank, NUM_OF_DIMS, cords);
	
	
	//printf("I'm a %d my cords: (%d, %d)\n", rank, cords[0], cords[1]);
	MPI_Comm MPI_COMM_ROW;
	MPI_Comm MPI_COMM_COL;

	MPI_Cart_sub(MPI_COMM_GRID, (int[]){0, 1}, &MPI_COMM_COL);
	MPI_Cart_coords(MPI_COMM_COL, rank, NUM_OF_DIMS, cords);
	printf("I'm a %d my cords: (%d, %d)\n", rank, cords[0], cords[1]);
	MPI_Finalize();
	return 0;
}
