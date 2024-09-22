#include <cstddef>
#include <cstdlib>
#include <iostream>

#include <math.h>
#include <mpi.h>
#include <pthread.h>

#define LISTS_COUNT  5
#define WEIGHT_COEFF 5000

#define MIN_COUNT_TASKS_TO_SHARE 20
#define TASK_PER_PROC            2400

#define TAG_REQUEST  0
#define TAG_RESPONSE 1

double RES_PER_ITER   = 0;
double GLOBAL_RES_SIN = 0;

int rank, size;

int* tasks;
int tasksRemain;
int tasksExecuted;

pthread_mutex_t mutexTasksList;
pthread_mutex_t mutexTasksRemain;

pthread_t recvThread;

void calculateTasks() {
	pthread_mutex_lock(&mutexTasksRemain);

	for (int i = 0; tasksRemain; ++i, tasksRemain--) {
		pthread_mutex_unlock(&mutexTasksRemain);


		pthread_mutex_lock(&mutexTasksList);
		int taskWeight = tasks[i];
		pthread_mutex_unlock(&mutexTasksList);

		for (int j = 0; j < taskWeight; ++j) {
			RES_PER_ITER += sin(j);
		}
		tasksExecuted++;

		pthread_mutex_lock(&mutexTasksRemain);
	}
	pthread_mutex_unlock(&mutexTasksRemain);
}

void initTasksWeight(void) {
	pthread_mutex_lock(&mutexTasksList);
	for (int i = 0; i < TASK_PER_PROC; ++i) {
		tasks[i] = abs(50 - i % 100) * abs(rank - (TASK_PER_PROC % size)) * WEIGHT_COEFF;
	}
	pthread_mutex_unlock(&mutexTasksList);
}

void* receiverRoutine(void* args) {
	int tasksToSend;
	int requesterRank;

	while (true) {
		//Get proc rank, that requests tasks
		MPI_Recv(&requesterRank, 1, MPI_INT, MPI_ANY_SOURCE, TAG_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		//Finish condition
		if (requesterRank == rank) {
			break;
		}

		pthread_mutex_lock(&mutexTasksRemain);
		if (tasksRemain >= MIN_COUNT_TASKS_TO_SHARE) {
			tasksToSend = tasksRemain / 2;
			tasksRemain -= tasksToSend;

			//Send tasks count
			MPI_Send(&tasksToSend, 1, MPI_INT, requesterRank, TAG_RESPONSE, MPI_COMM_WORLD);

			pthread_mutex_lock(&mutexTasksList);
			//Send tasks
			MPI_Send(tasks + tasksExecuted + tasksRemain - 1, tasksToSend, MPI_INT, requesterRank, TAG_RESPONSE, MPI_COMM_WORLD);
			
			pthread_mutex_unlock(&mutexTasksRemain);
			pthread_mutex_unlock(&mutexTasksList);

		} else {
			tasksToSend = 0;
			MPI_Send(&tasksToSend, 1, MPI_INT, requesterRank, TAG_RESPONSE, MPI_COMM_WORLD);
		}
	}

	return NULL;
}

void* workerRoutine(void* args) {
	tasks = new int[TASK_PER_PROC];

	double startTime;
	double minTime, maxTime;

	for (int iterCount = 0; iterCount < LISTS_COUNT; ++iterCount) {
		initTasksWeight();

		pthread_mutex_lock(&mutexTasksRemain);
		tasksRemain = TASK_PER_PROC;
		pthread_mutex_unlock(&mutexTasksRemain);

		tasksExecuted       = 0;
		int additionalTasks = 0;

		startTime = MPI_Wtime();
	
		//Work on his own tasks
		calculateTasks();

		//Do someone else work
		for (int currProc = 0; currProc < size; ++ currProc) {
			if (currProc == rank) {
				continue;
			}

			//Peek message about readiness to accept a work
			MPI_Send(&rank, 1, MPI_INT, currProc, TAG_REQUEST, MPI_COMM_WORLD);

			//Get tasks count
			MPI_Recv(&additionalTasks, 1, MPI_INT, currProc, TAG_RESPONSE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			//If there are tasks, get them and start execute
			if (additionalTasks > 0) {
				MPI_Recv(tasks, additionalTasks, MPI_INT, currProc, TAG_RESPONSE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				pthread_mutex_lock(&mutexTasksRemain);
				tasksRemain = additionalTasks;
				pthread_mutex_unlock(&mutexTasksRemain);

				calculateTasks();
			}
		}

		double endTime = MPI_Wtime();
		double resTime = endTime - startTime;

		MPI_Allreduce(&resTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&resTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		if (rank== 0) {
			std::cout << "=================================================" << std::endl;

			std::cout << "Iteration numer: " << iterCount << std::endl;
			std::cout << "Disbalance time: " << maxTime - minTime << std::endl;
			std::cout << "Disbalance percentage: " << (maxTime - minTime) / maxTime * 100 << std::endl;
			std::cout << "----------------------------------------------"  << std::endl;
		}

		for (int currProc = 0; currProc < size; currProc++) {
			if (rank == currProc) {
				std::cout << "\t\tCurrent proc is: " << rank << std::endl;
				std::cout << "Amount of executed tasks: " << tasksExecuted << std::endl;
				std::cout << "Result of calculating is: " << RES_PER_ITER << std::endl;
				std::cout << "Time per iteration: " << resTime << " seconds" << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		
			
	}
	//Tell receiver about finish
	MPI_Send(&rank, 1, MPI_INT, rank, TAG_REQUEST, MPI_COMM_WORLD);

	MPI_Allreduce(&RES_PER_ITER, &GLOBAL_RES_SIN, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	delete[] tasks; 
	return NULL;
}

void startThreads(void) {
	pthread_mutex_init(&mutexTasksList, NULL);
	pthread_mutex_init(&mutexTasksRemain, NULL);

	if (pthread_create(&recvThread, NULL, receiverRoutine, NULL) != 0) {
		MPI_Finalize();
		perror("Can't create thread");
		exit(EXIT_FAILURE);
	}

	workerRoutine(NULL);

	pthread_mutex_destroy(&mutexTasksList);
	pthread_mutex_destroy(&mutexTasksRemain);
}

int main(int argc, char** argv) {
	int reqLevel  = MPI_THREAD_MULTIPLE;
	int provLevel;

	MPI_Init_thread(&argc, &argv, reqLevel, &provLevel);
	if (provLevel != reqLevel) {
		MPI_Finalize();
		perror("Can't load reqLevel");
		return 0;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double startTime = MPI_Wtime();
	startThreads();
	double endTime = MPI_Wtime();

	if (rank == 0) {
		std::cout << "Time for all lists: " << endTime - startTime << std::endl; 
	}

	MPI_Finalize();
	return 0;
}
