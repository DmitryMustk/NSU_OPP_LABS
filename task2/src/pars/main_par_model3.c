#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void fillMatrix(double* matrix, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		srand(i);
		for (int j = 0; j < vecSize; ++j) {
			if (i == j) {
				matrix[i * vecSize + j] = rand() % (vecSize / 3) + vecSize;
				continue;
			}
			matrix[i * vecSize + j] = vecSize + 1;
		}
	}
}

void fillVecU(double* vec, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		vec[i] = sin((2 * M_PI * i) / vecSize);
	}
}


void printMatrix(double* matrix, int rows, int cols, int rank) {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			printf("%.1f ", matrix[i * cols + j]);
		}
		printf("\n");
	}
	printf("-- from %d process\n", rank);
}


void fillVecB(double* vec, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		vec[i] = vecSize + 1;
	}
}

void mulMatrixOnVec(double* matrix, double* vec, double* res, int rows, int vecSize) {
	for (int i = 0; i < rows; ++i) {
		res[i] = 0;
		for (int j = 0; j < vecSize; ++j) {
			res[i] += matrix[i * vecSize + j] * vec[j];
		}
	}
}

void vectorSub(double* vec1, double* vec2, double* res, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		res[i] = vec1[i] - vec2[i];
	}
}

void vectorSum(double* vec1, double* vec2, double* res, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		res[i] = vec1[i] + vec2[i];
	}
}

void copyVector(double* fromVec, double* toVec, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		toVec[i] = fromVec[i];
	}
}

double vectorDotProduct(double* vec1, double* vec2, int vecSize) {
	double res = 0;
	for (int i = 0; i < vecSize; ++i) {
		res += vec1[i] * vec2[i];
	}
	return res;
}

void mulVecOnScalar(double* vec, double* res, double scalar, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		res[i] = vec[i] * scalar;
	}
}

double getVecSquaredNorm(double* vec, int vecSize) {
	double sqNorm = 0;
	for (int i = 0; i < vecSize; ++i) {
		sqNorm += vec[i] * vec[i];
	}
	return sqNorm;
}

void printVec(double* vec, int vecSize, int rank) {
	printf("(");
	for (int i = 0; i < vecSize; ++i) {
		printf("%f, ", vec[i]);
	}
	printf(") -- from %d process\n", rank);
}


//for Scatterv int* sendcount
void getElemsOnProc(int* elemsOnProc, int vecSize, int size) {
	int baseRowOnProc = vecSize / size;
	int restRows = vecSize % size;

	for (int i = 0; i < size; ++i) {
		elemsOnProc[i] = baseRowOnProc * vecSize;

		if (restRows > 0) {
			elemsOnProc[i] += vecSize;
			restRows--;
		}
	}
}

//for Scatterv int* displs
void getElemsOffset(int* elemsOffset, int* elemsOnProc, int size) {
	int offset = 0;

	for (int i = 0; i < size; ++i) {
		elemsOffset[i] = offset;
		offset += elemsOnProc[i];
	}
}

void getRowsOnProc(int* rowsOnProc, int vecSize, int size) {
	int baseRowOnProc = vecSize / size;
	int restRows = vecSize % size;

	for (int i = 0; i < size; ++i) {
		rowsOnProc[i] = baseRowOnProc;

		if (restRows > 0) {
			rowsOnProc[i]++;
			restRows--;
		}
	} 
}

void getRowsOffset(int* rowsOffset, int* rowsOnProc, int size) {
	int offset = 0;

	for (int i = 0; i < size; ++i) {
		rowsOffset[i] = offset;
		offset += rowsOnProc[i];
	}
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int vecSize = 10000;	

	double* matrix  = malloc(sizeof(double) * vecSize * vecSize);
	double* vecB    = malloc(sizeof(double) * vecSize);
	double* vecU    = malloc(sizeof(double) * vecSize);
	
	double  epsSq   = 0.000001 * 0.000001;

	double* vecX    = calloc(vecSize, sizeof(double));
	double* vecR    = malloc(sizeof(double) * vecSize);
	double* vecZ    = malloc(sizeof(double) * vecSize);

	double* tmp1Vec = malloc(sizeof(double) * vecSize);
	double* tmp2Vec = malloc(sizeof(double) * vecSize);

	double alpha;
	double betta;

	//
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int* rowsOnProc   = malloc(sizeof(int) * size);
	int* rowsOffset   = malloc(sizeof(int) * size);
	int* elemsOnProc  = malloc(sizeof(int) * size);
	int* elemsOffset  = malloc(sizeof(int) * size);

	getRowsOnProc(rowsOnProc, vecSize, size);
	getRowsOffset(rowsOffset, rowsOnProc, size);
	getElemsOnProc(elemsOnProc, vecSize, size);
	getElemsOffset(elemsOffset, elemsOnProc, size);

	double* portionMatrix = malloc(sizeof(double) * elemsOnProc[rank]);


	if (rank == 0) {
		fillMatrix(matrix, vecSize);
		fillVecU(vecU, vecSize);
		mulMatrixOnVec(matrix, vecU, vecB, vecSize, vecSize);
	}

	double start = 0;
	if (rank == 0) {
		start = MPI_Wtime();
	}

	MPI_Bcast(vecB, vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(matrix, elemsOnProc, elemsOffset, MPI_DOUBLE, portionMatrix, elemsOnProc[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double vecBSqNorm;
	double localVecBSqNorm = getVecSquaredNorm(vecB + rowsOffset[rank], rowsOnProc[rank]);
	MPI_Allreduce(&localVecBSqNorm, &vecBSqNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   //r^0
	mulMatrixOnVec(portionMatrix, vecX, tmp1Vec, rowsOnProc[rank], vecSize);
	vectorSub(vecB + rowsOffset[rank], tmp1Vec, tmp2Vec, rowsOnProc[rank]);
	MPI_Allgatherv(tmp2Vec, rowsOnProc[rank], MPI_DOUBLE, vecR, rowsOnProc, rowsOffset, MPI_DOUBLE, MPI_COMM_WORLD);
	//printVec(vecR, vecSize, rank);
	
	//z^0
	copyVector(vecR, vecZ, vecSize);

	int matches = 0;
	int i;
	double rSqScalar;
	double Azz;
	double vecRSqNorm;
	for (i = 0; i < 10000; ++i) {
		double localRSqScalar = vectorDotProduct(vecR + rowsOffset[rank], vecR + rowsOffset[rank], rowsOnProc[rank]);
		MPI_Allreduce(&localRSqScalar, &rSqScalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		//alpha
		mulMatrixOnVec(portionMatrix, vecZ, tmp1Vec, rowsOnProc[rank], vecSize); //tmp1 -> A*z
		double localAzz = vectorDotProduct(tmp1Vec, vecZ + rowsOffset[rank], rowsOnProc[rank]);
		MPI_Allreduce(&localAzz, &Azz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		alpha = rSqScalar / Azz;
		
		//x^n+1
		mulVecOnScalar(vecZ + rowsOffset[rank], tmp2Vec, alpha, rowsOnProc[rank]); //alpha * z^n
		vectorSum(vecX + rowsOffset[rank], tmp2Vec, tmp2Vec, rowsOnProc[rank]);
		MPI_Allgatherv(tmp2Vec, rowsOnProc[rank], MPI_DOUBLE, vecX, rowsOnProc, rowsOffset, MPI_DOUBLE, MPI_COMM_WORLD); //

		//r^n+1
		mulVecOnScalar(tmp1Vec, tmp2Vec, alpha, rowsOnProc[rank]);
		vectorSub(vecR + rowsOffset[rank], tmp2Vec, tmp1Vec, rowsOnProc[rank]); //tmp1 -> portion of R
		MPI_Allgatherv(tmp1Vec, rowsOnProc[rank], MPI_DOUBLE, vecR, rowsOnProc, rowsOffset, MPI_DOUBLE, MPI_COMM_WORLD); //

		double localVecRSqNorm = getVecSquaredNorm(tmp1Vec, rowsOnProc[rank]);
		MPI_Allreduce(&localVecRSqNorm, &vecRSqNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if (vecRSqNorm / vecBSqNorm < epsSq) { //-------------------------exit condition
			matches++;
			if (matches == 3) {
				break;
			}
		} else {
			matches = 0;
		}
		
		//betta
		double localBetta = vectorDotProduct(tmp1Vec, tmp1Vec, rowsOnProc[rank]);
		MPI_Allreduce(&localBetta, &betta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		betta  /= rSqScalar;
		//z^n+1
		mulVecOnScalar(vecZ + rowsOffset[rank], tmp1Vec, betta, rowsOnProc[rank]);
		vectorSum(vecR + rowsOffset[rank], tmp1Vec, tmp2Vec, rowsOnProc[rank]);
		MPI_Allgatherv(tmp2Vec, rowsOnProc[rank], MPI_DOUBLE, vecZ, rowsOnProc, rowsOffset, MPI_DOUBLE, MPI_COMM_WORLD); //
	}
	if (rank == 0) {
		printf("done on %d iter. Time taken: %f\n", i, MPI_Wtime() - start);
		printVec(vecX, vecSize, 0);
	}
	free(matrix);
	free(vecB);
	free(vecU);
	free(vecX);
	free(tmp1Vec);
	free(tmp2Vec);
	free(vecR);
	free(vecZ);

	MPI_Finalize();
	return 0;
}
