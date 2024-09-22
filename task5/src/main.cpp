#include <cstdlib>
#include <mpi.h>

#include <iostream>
#include <ostream>
#include <vector>

void generateGlider(bool* matrix, int colsCount, int rowsCount) {
	std::fill(matrix, matrix + colsCount * rowsCount, false);

	matrix[0 * colsCount + 1] = true;
	
	matrix[1 * colsCount + 2] = true;
	
	matrix[2 * colsCount + 0] = true;
	matrix[2 * colsCount + 1] = true;
	matrix[2 * colsCount + 2] = true;
}

int countNbours(bool *oldData, int columnsAmount, int i, int j) {
  int neighboirsAmount =
      (oldData[i * columnsAmount + (j + 1) % columnsAmount]) +
      (oldData[i * columnsAmount +
               (j + columnsAmount - 1) % columnsAmount]) +
      (oldData[(i + 1) * columnsAmount + (j + 1) % columnsAmount]) +
      (oldData[(i + 1) * columnsAmount +
               (j + columnsAmount - 1) % columnsAmount]) +
      (oldData[(i - 1) * columnsAmount + (j + 1) % columnsAmount]) +
      (oldData[(i - 1) * columnsAmount +
               (j + columnsAmount - 1) % columnsAmount]) +
      (oldData[(i + 1) * columnsAmount + j]) +
      (oldData[(i - 1) * columnsAmount + j]);

  return neighboirsAmount;
}

void computeNextGeneration(bool* oldMatrix, bool* nextMatrix, int rowsCount, int colsCount) {
	for (int i = 1; i < rowsCount - 1; ++i) {
		for (int j = 0; j < colsCount; ++j) {
			int state = oldMatrix[i * colsCount + j];

			int nboursCount = countNbours(oldMatrix, colsCount, i, j);

			if (state == 0 && nboursCount == 3) {
				nextMatrix[i * colsCount + j] = 1;
			}
			else if (state == 1 && (nboursCount < 2 || nboursCount > 3)) {
				nextMatrix[i * colsCount + j] = 0;
			}
			else {
				nextMatrix[i * colsCount + j] = oldMatrix[i * colsCount + j] = state;
			}
		}
	}
}

bool isStop(bool* stopMatrix, int rowsCount, int colsCount) {
	for (int i = 0; i < colsCount; ++i) {
		bool stop = true;
		for (int j = 0; j < rowsCount; ++j) {
			stop &= stopMatrix[j * colsCount + i];
		}
		if (stop) {
			return true;
		}
	}

	return false;
} 

bool matEqual(bool* matrix1, bool* matrix2, int realStart, int realEnd) {
	for (int i = realStart; i < realEnd; ++i) {
		if (matrix1[i] != matrix2[i]) {
			return false;
		}
	}
	return true;
}

void calcStopVector(std::vector<bool*> worldHistory, bool* stopVector, bool* extendedPortionMatrix, int rowsCount, int colsCount) {

	int vecSize = worldHistory.size() - 1;
	auto it = worldHistory.begin();
	for (int i = 0; i < vecSize; ++i) {
		stopVector[i] = matEqual(*it, extendedPortionMatrix, colsCount, colsCount * (rowsCount + 1));
		it++;
	}
} 

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

void startLife(bool* startMatrix, int rowsCount, int colsCount, int rank, int size) {

	int* rowsOnProc  = new int[size];
	int* rowsOffset  = new int[size];
	int* elemsOnProc = new int[size];
	int* elemsOffset = new int[size];

	getRowsOnProc(rowsOnProc, rowsCount, size);
	getRowsOffset(rowsOffset, rowsOnProc, size);
	getElemsOnProc(elemsOnProc, rowsCount, size);
	getElemsOffset(elemsOffset, elemsOnProc, size);

	bool* extendedPortionMatrix = new bool[elemsOnProc[rank] + colsCount * 2];
	bool* basePortionMatrix     = extendedPortionMatrix + colsCount;

	MPI_Scatterv(startMatrix, elemsOnProc, elemsOffset, MPI_C_BOOL, basePortionMatrix, elemsOnProc[rank], MPI_C_BOOL, 0, MPI_COMM_WORLD);

	int prevRank = (rank + size - 1) % size;
	int nextRank = (rank + 1) % size;

	std::vector<bool*> worldHistory; 

	int  iter = 0;
	bool stop = false;

	while (!stop) {
		bool* nextExtendedPortionMatrix = new bool[elemsOnProc[rank] + 2 * colsCount];
		bool* nextBasePortionMatrix = nextExtendedPortionMatrix + colsCount;
		
		worldHistory.push_back(extendedPortionMatrix);
		iter++;
		
		MPI_Request sendFirstLineReq;
		MPI_Request sendLastLineReq;
		// 1 - initiation of sending first line to the prev proc
		MPI_Isend(basePortionMatrix, colsCount, MPI_C_BOOL, prevRank, 1, MPI_COMM_WORLD, &sendFirstLineReq);
		// 2 - initiation of sending last line to the next proc
		MPI_Isend(basePortionMatrix + elemsOnProc[rank] - colsCount, colsCount, MPI_C_BOOL, nextRank, 2, MPI_COMM_WORLD, &sendLastLineReq);


		MPI_Request recvLastLineReq;
		MPI_Request recvFirstLineReq;
		// 3 - initiation of receiving last line from the prev proc
		MPI_Irecv(extendedPortionMatrix, colsCount, MPI_C_BOOL, prevRank, 2, MPI_COMM_WORLD, &recvLastLineReq);
		// 4 - initiation of receiveng first line from the next proc
		MPI_Irecv(basePortionMatrix + elemsOnProc[rank], colsCount, MPI_C_BOOL, nextRank, 1, MPI_COMM_WORLD, &recvFirstLineReq);


		// 5 - count vector of stop flags
		MPI_Request flagsReq;
		bool* stopVector;
		bool* stopMatrix;

		int vecSize = worldHistory.size() - 1;
		if (vecSize > 1) {
			stopVector = new bool[vecSize];
			calcStopVector(worldHistory, stopVector, extendedPortionMatrix, rowsOnProc[rank], colsCount);	

			stopMatrix = new bool[vecSize * size];

			// 6 - initiation gathering stop vectors
			MPI_Iallgather(stopVector, vecSize, MPI_C_BOOL, stopMatrix, vecSize, MPI_C_BOOL, MPI_COMM_WORLD, &flagsReq);
		}
		
		// 7 - compute base mat
		computeNextGeneration(basePortionMatrix, nextBasePortionMatrix, rowsOnProc[rank], colsCount);


		MPI_Status status;

		// 8 - wait sending first line to prev proc
		MPI_Wait(&sendFirstLineReq, &status);
		
		// 9 - wait end of receiveng last line from the next proc
		MPI_Wait(&recvLastLineReq, &status);
		
		// 10 - count states of the first line (need to give first and get last)
		computeNextGeneration(extendedPortionMatrix, nextExtendedPortionMatrix, 3, colsCount);


		// 11 - wait sending last line to next proc
		MPI_Wait(&sendLastLineReq, &status);

		// 12 - wait receiveng first line
		MPI_Wait(&recvFirstLineReq, &status);
		
		// 13 - count stages of the last line (need to give last and get first)
		computeNextGeneration(
				basePortionMatrix + (rowsOnProc[rank] - 2) * colsCount,
				nextBasePortionMatrix + (rowsOnProc[rank] - 2), 3, colsCount);

		if (vecSize > 1) {
			// 14 - wait receiveng other stop vectors
			MPI_Wait(&flagsReq, &status);

			// 15 - compare stop vectors
			stop = isStop(stopMatrix, size, vecSize);


			delete[] stopVector;
			delete[] stopMatrix;
		}

		if (stop) {
			break;
		}

		extendedPortionMatrix = nextExtendedPortionMatrix;
		basePortionMatrix     = nextBasePortionMatrix;
	}
	
	if (rank == 0) {
		std::cout << "Repeat after: " << iter - 1 << "iterations" << std::endl;
	}

	for (auto mat : worldHistory) {
		delete[] mat;
	}
	
	delete[] rowsOnProc;
	delete[] rowsOffset;
	delete[] elemsOnProc;
	delete[] elemsOffset;
	
}

int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "Need to pass <rowsCount> <colsCount>" << std::endl;
		return 0;
	}

	int rowsCount = std::atoi(argv[1]);
	int colsCount = std::atoi(argv[2]);

	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	bool* startMatrix = nullptr;
	if (rank == 0) {
		startMatrix = new bool[rowsCount * colsCount];
		generateGlider(startMatrix, colsCount, rowsCount);
	}

	double startTime = MPI_Wtime();
	startLife(startMatrix, rowsCount, colsCount, rank, size);
	double endTime = MPI_Wtime();

	if (rank == 0) {
		std::cout << "Time taken: " << endTime - startTime << std::endl;
		delete[] startMatrix;
	}

	MPI_Finalize();
	return 0;
}
