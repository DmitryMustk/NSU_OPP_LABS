#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void fillMatrix(double* matrix, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		for (int j = 0; j < vecSize; ++j) {
			if (i == j) {
				matrix[i * vecSize + j] = 2.0;
				continue;
			}
			matrix[i * vecSize + j] = 1.0;
		}
	}
}

void printMatrix(double* matrix, int raws, int cols) {
	for (int i = 0; i < raws; ++i) {
		for (int j = 0; j < cols; ++j) {
			printf("%.1f ", matrix[i * cols + j]);
		}
		printf("\n");
	}
}


void fillVecB(double* vec, int vecSize) {
	for (int i = 0; i < vecSize; ++i) {
		vec[i] = vecSize + 1;
	}
}

void mulMatrixOnVec(double* matrix, double* vec, double* res, int vecSize) {
	#pragma omp for
	for (int i = 0; i < vecSize; ++i) {
		res[i] = 0;
		for (int j = 0; j < vecSize; ++j) {
			res[i] += matrix[i * vecSize + j] * vec[j];
		}
	}
}

void vectorSub(double* vec1, double* vec2, double* res, int vecSize) {
	#pragma omp for
	for (int i = 0; i < vecSize; ++i) {
		res[i] = vec1[i] - vec2[i];
	}
}

void vectorSum(double* vec1, double* vec2, double* res, int vecSize) {
	#pragma omp for
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
	#pragma omp parallel for reduction(+:res)
	for (int i = 0; i < vecSize; ++i) {
		res += vec1[i] * vec2[i];
	}
	return res;
}

void mulVecOnScalar(double* vec, double* res, double scalar, int vecSize) {
	#pragma omp for
	for (int i = 0; i < vecSize; ++i) {
		res[i] = vec[i] * scalar;
	}
}

double getVecSquaredNorm(double* vec, int vecSize) {
	double sqNorm = 0;
	#pragma omp parallel for reduction(+:sqNorm)
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

int main(void) {
	int vecSize = 1000;	

	double* matrix = malloc(sizeof(double) * vecSize * vecSize);
	double* vecB = malloc(sizeof(double) * vecSize);
	double epsSquared = 0.000001 * 0.000001;

	double* vecX = calloc(vecSize, sizeof(double));
	double* tmp1Vec = malloc(sizeof(double) * vecSize);
	double* tmp2Vec = malloc(sizeof(double) * vecSize);
	double* vecR = malloc(sizeof(double) * vecSize);
	double* vecZ = malloc(sizeof(double) * vecSize);

	double alpha;
	double betta;

	fillMatrix(matrix, vecSize);
	fillVecB(vecB, vecSize);

	int matches = 0;

	
	#pragma omp parallel private(matches)
	{
		
		double vecBSqNorm = getVecSquaredNorm(vecB, vecSize);
		int rank = omp_get_thread_num();
		//r^0
		mulMatrixOnVec(matrix, vecX, tmp1Vec, vecSize);
		vectorSub(vecB, tmp1Vec, vecR, vecSize);
		//printf("r_0 norm: %f", getVecSquaredNorm(vecR, vecSize));

		//z^0
		copyVector(vecR, vecZ, vecSize);

		int i;
		for (i = 0; i < 100; ++i) {
			//printVec(vecR, vecSize, 0);
		
			double rSqScalar = vectorDotProduct(vecR, vecR, vecSize);
			//alpha
			mulMatrixOnVec(matrix, vecZ, tmp1Vec, vecSize);
			alpha = rSqScalar / vectorDotProduct(tmp1Vec, vecZ, vecSize);
			//x^n+1
			mulVecOnScalar(vecZ, tmp2Vec, alpha, vecSize); //alpha * z^n
			vectorSum(vecX, tmp2Vec, vecX, vecSize);
			//printVec(vecX, vecSize, 0);
			//r^n+1
			mulVecOnScalar(tmp1Vec, tmp2Vec, alpha, vecSize);
			vectorSub(vecR, tmp2Vec, vecR, vecSize);
		
			//printVec(vecR, vecSize, 0);
			double vecRSqNorm = getVecSquaredNorm(vecR, vecSize);
			if (vecRSqNorm / vecBSqNorm < epsSquared) {
				//printf("iter count: %d\n", i);
				matches++;
				if (matches == 1) {
					break;
				}
				//printf("r norm: %f, b norm: %f\n",vecRSqNorm, vecBSqNorm);
			} else {
				matches = 0;
			}
		
			//betta
			betta = vectorDotProduct(vecR, vecR, vecSize) / rSqScalar;
			//printf("||betta:%f", betta);
			//z^n+1
			mulVecOnScalar(vecZ, tmp1Vec, betta, vecSize); // betta * z^n
			vectorSum(vecR, tmp1Vec, vecZ, vecSize);
			//printVec(vecZ, vecSize, 0);

		}
		if (rank == 0) {
			printf("done on %d iter\n", i);
			printVec(vecX, vecSize, 0);
		}
	}
	free(matrix);
	free(vecB);
	free(vecX);
	free(tmp1Vec);
	free(tmp2Vec);
	free(vecR);
	free(vecZ);

	return 0;
}
