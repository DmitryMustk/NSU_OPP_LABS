#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>

void printSquareMatrix(double* matrix, int matrixSize) {
	for (int i = 0; i < matrixSize; ++i) {
		for (int j = 0; j < matrixSize; ++j) {
			printf("%.2f ", matrix[i * matrixSize + j]);
		}
		printf("\n");
	}
}

void printMatrix(double* matrix, int raws, int cols) {
	for (int i = 0; i < raws; ++i) {
		for (int j = 0; j < cols; ++j) {
			printf("%.0f ", matrix[i * cols + j]);
		}
		printf("\n");
	}
}

void fillMatrix(double* matrix, int matrixSize) {
	for (int i = 0; i < matrixSize; ++i) {
		for (int j = 0; j < matrixSize; ++j) {
			matrix[i * matrixSize + j] = i;
		}
	}
}

int* getElementsOnProc(int vecSize, int size) {
	int rowOnProc = vecSize / size;
	int restRows = vecSize % size;

	int* elemsCountArr = malloc(sizeof(int) * size);

	for (int i = 0; i < size; ++i) {
		elemsCountArr[i] = rowOnProc * vecSize;

		if (restRows > 0) {
			elemsCountArr[i] += vecSize;
			restRows--;
		}
	} 
	return elemsCountArr;
}

//for Scatterv int* sendcounts
int* getMatrixSendCounts(int vecSize, int size) {
	int* arr = malloc(sizeof(int) * size);

	int baseRawCount = vecSize / size;
	int restRawCount = vecSize % size;

	for (int i = 0; i < size; ++i) {
		arr[i] = baseRawCount * vecSize;

		if (restRawCount > 0) {
			arr[i] += vecSize;
			restRawCount--;
		}
	}
	return arr;
}

//for Scatterv int* displs
int* getMatrixDispls(int* matrixSendCounts, int size) {
	int* arr = malloc(sizeof(int) * size);

	int offset = 0;
	for (int i = 0; i < size; ++i) {
		arr[i] = offset;
		offset += matrixSendCounts[i];
	}

	return arr;
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int vecSize = 10;
	int* sendCountsArr = getMatrixSendCounts(vecSize, size);
	int* displsArr = getMatrixDispls(sendCountsArr, size);

	double* matrix = malloc(sizeof(double) * vecSize * vecSize);
	double* portionMatrix = malloc(sizeof(double) * sendCountsArr[rank]);
	//double* portionMatrix = malloc(sizeof(double) * vecSize / size * vecSize);

	if (rank == 0) {
		fillMatrix(matrix, vecSize);
	}	
	
	//MPI_Scatter(matrix, vecSize * vecSize, MPI_DOUBLE, portionMatrix, vecSize / size * vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(matrix, sendCountsArr, displsArr, MPI_DOUBLE, portionMatrix, sendCountsArr[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 2) {
		printMatrix(portionMatrix, sendCountsArr[rank] / vecSize, vecSize);
		//printMatrix(matrix, vecSize, vecSize);
	}
	
	MPI_Finalize();
}
