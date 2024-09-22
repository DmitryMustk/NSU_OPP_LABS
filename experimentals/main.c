#include <mpi.h>
#include <stdio.h>

int main(void) {
	printf("%ld", sizeof(double));
	printf("%ld", MPI_DOUBLE);
	return 0;
}
