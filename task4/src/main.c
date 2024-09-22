#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ERROR -1

#define NUM_OF_DIMS 2

#define ROW_CORD 0
#define COL_CORD 1

void printMatrix(double* matrix, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			printf("%f ", matrix[i * cols +j]);
		}
		printf("\n");
	}
} 

int main(int argc, char** argv) {
	if (argc < 3) {
		printf("Need to pass 2 args");
		return ERROR;
	}

	int p1 = atoi(argv[1]);
	int p2 = atoi(argv[2]); 

	int dims[] = {p1, p2};
	int periods[] = {0, 0};

	int cords[NUM_OF_DIMS];
	
	int n1 = 6168;
	int n2 = 3000;
	int n3 = 5184;

	int rowsInSlice = n1 / p1;
	int colsInSlice = n3 / p2;

	double* matrixA;
	double* matrixB;
	double* matrixRes;

	
	double* portionMatA   = malloc(sizeof(double) * rowsInSlice * n2);
	double* portionMatB   = malloc(sizeof(double) * n2 * colsInSlice);
	double* portionMatRes = malloc(sizeof(double) * rowsInSlice * colsInSlice);

	double startTime;
	
	MPI_Init(NULL, NULL);
	
	int rank, size;
	int colSize, rowSize;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm MPI_COMM_GRID;
	MPI_Cart_create(MPI_COMM_WORLD, NUM_OF_DIMS, dims, periods, 0, &MPI_COMM_GRID);
	

	MPI_Cart_coords(MPI_COMM_GRID, rank, NUM_OF_DIMS, cords);

	MPI_Comm MPI_COMM_ROWS;
	MPI_Comm MPI_COMM_COLS;

	MPI_Cart_sub(MPI_COMM_GRID, (int[]){0, 1}, &MPI_COMM_ROWS);
	MPI_Cart_sub(MPI_COMM_GRID, (int[]){1, 0}, &MPI_COMM_COLS);

	MPI_Comm_size(MPI_COMM_ROWS, &rowSize);
	MPI_Comm_size(MPI_COMM_ROWS, &colSize);

	if (cords[ROW_CORD] == 0 && cords[COL_CORD] == 0) {
		printf("grid size = %dx%d\n", rowSize, colSize);
	}

	//INITIALIZE
	if (cords[ROW_CORD] == 0 && cords[COL_CORD] == 0) {
		matrixA   = malloc(sizeof(double) * n1 * n2);
		matrixB   = malloc(sizeof(double) * n2 * n3);
		matrixRes = malloc(sizeof(double) * n1 * n3);
		for (int i = 0; i < n1 * n2; ++i) {
			matrixA[i] = 1.0;
		}
		for (int i = 0; i < n2 * n3; ++i) {
			matrixB[i] = 2.0;
		}
		startTime = MPI_Wtime();
	}
	//MATRIX A DISTIBUTION
		//(0, 0) -> (1, 0), (2, 0), ...
	if (cords[COL_CORD] == 0) {
		MPI_Scatter(matrixA, rowsInSlice * n2, MPI_DOUBLE, portionMatA, rowsInSlice * n2, MPI_DOUBLE, 0, MPI_COMM_COLS);
	}
		//(0, 0) -> (0, 1), (0, 2), ...
		//(1, 0) -> (1, 1), (1, 2), ...
		//.............................
	MPI_Bcast(portionMatA, rowsInSlice * n2, MPI_DOUBLE, 0, MPI_COMM_ROWS);
	
	//MATRIX B DISTIBUTION
	MPI_Datatype MPI_MAT_B_COL;
	MPI_Type_vector(n2, colsInSlice, n3, MPI_DOUBLE, &MPI_MAT_B_COL);
	MPI_Type_commit(&MPI_MAT_B_COL);
	
		//(0, 0) -> (0, 1), (0, 2), ...
	if (cords[ROW_CORD] == 0 && cords[COL_CORD] == 0) {
		//(0,0) gives matB to himself
		for (int i = 0; i < n2; ++i) {
			for (int j = 0; j < colsInSlice; ++j) {
				portionMatB[i * colsInSlice + j] = matrixB[i * n3 + j];
			}
		}
		//(0,0) gives matB to others
		for (int i = 1; i < rowSize; ++i) {
			MPI_Send(matrixB + i * colsInSlice, 1, MPI_MAT_B_COL, i, 0, MPI_COMM_ROWS);
		}
	}
	else if (cords[ROW_CORD] == 0) {
		MPI_Recv(portionMatB, colsInSlice * n2, MPI_DOUBLE, 0, 0, MPI_COMM_ROWS, MPI_STATUS_IGNORE);
	}
		//(0, 0) -> (1, 0), (2, 0), ...
		//(0, 1) -> (1, 1), (2, 1), ...
		//.............................
	MPI_Bcast(portionMatB, n2 * colsInSlice, MPI_DOUBLE, 0, MPI_COMM_COLS);

	
	//PORTION_A X POTION_B -> PORTION_C
	for (int i = 0; i < rowsInSlice; ++i) {
		for (int j = 0; j < colsInSlice; ++j) {
			portionMatRes[i * colsInSlice + j] = 0.0;
			for (int k = 0; k < n2; ++k) {
				portionMatRes[i * colsInSlice + j] +=
					portionMatA[i * n2 + k] * portionMatB[k * colsInSlice + j];
			}
		}
	}

	printf("portion C on (%d, %d) process\n", cords[0], cords[1]);
	printMatrix(portionMatRes, rowsInSlice, colsInSlice);
	
	//portionMatRes[i] -> matrixRes
	MPI_Datatype MPI_PORTION_RES;
	MPI_Type_vector(rowsInSlice, colsInSlice, n3, MPI_DOUBLE, &MPI_PORTION_RES);
	MPI_Type_commit(&MPI_PORTION_RES);
	
	if (cords[ROW_CORD] == 0 && cords[COL_CORD] == 0) {
		//recieve matRes from others
		for (int i = 0; i < p1; ++i) {
			for (int j = 0; j < p2; ++j) {
				if (i == 0 && j == 0) {
					continue;
				}
				MPI_Recv(matrixRes + i * n3 * rowsInSlice + j * colsInSlice, 1, MPI_PORTION_RES, i * p2 + j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		//recieve matRes from himself
		for (int i = 0; i < rowsInSlice; ++i) {
			for (int j = 0; j < colsInSlice; ++j) {
				matrixRes[i * n3 + j] = portionMatRes[i * colsInSlice + j];
			}
		}
	}
	else {
		//send portionRes to (0, 0)
		MPI_Send(portionMatRes, rowsInSlice * colsInSlice, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	if (cords[ROW_CORD] == 0 && cords[COL_CORD] == 0) {
		printf("Time taken: %lf\n", MPI_Wtime() - startTime);
		printMatrix(matrixRes, n1, n3);
		free(matrixA);
		free(matrixB);
		free(matrixRes);
	}

	//free(portionMatA);
	//free(portionMatB);
	//free(portionMatRes);

	MPI_Type_free(&MPI_PORTION_RES);
	MPI_Type_free(&MPI_MAT_B_COL);
	
	MPI_Finalize();
	return 0;
}

