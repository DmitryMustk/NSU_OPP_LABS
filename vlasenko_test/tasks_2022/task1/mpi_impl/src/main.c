#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>

#define ERROR -1
#define TARGET_SUM 1000000
#define MAX_ITERS 10

int generate_number(void) {
    return rand() % 250000 + 100000;
}

void print_arr(int* arr, int size) {
    printf("|||");
    for(int i = 0; i < size; ++i) {
        printf("%d ", arr[i]);
    }
    printf("|||\n");
}

void null_arr(int* arr, int size) {
    for(int i = 0; i < size; ++i) {
        arr[i] = 0;
    }
}

int get_max_number(int* arr, int size) {
    int max = 0;
    for(int i = 0; i < size; ++i) {
        if(arr[i] > max) {
            max = arr[i];
        }
    }
    return max;
}

int get_min_number(int* arr, int size) {
    int min = INT_MAX;
    for(int i = 0; i < size; ++i) {
        if(arr[i] < min)
            min = arr[i];
    }
    return min;
}

int main(int argc, char *argv[]) {
    int rank, size;

    srand(time(NULL) % (rank + 1));

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) {
            printf("Usage: %s <number_of_processes>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }    

    int num_processes = atoi(argv[1]);

    if (num_processes > size) {
        if (rank == 0) {
            printf("Number of processes exceeds available MPI processes.\n");
        }
        MPI_Finalize();
        return ERROR;
    }

    int cur_number = 0;
    int cur_summ = 0;
    int cur_iter;

    int numbers[MAX_ITERS];
    int all_numbers[size * MAX_ITERS];
    int iters[size];
    int total_iters[size];
    null_arr(total_iters, size);

    int random_numbers[size];

    for(int lap = 1; lap <= 3; ++lap) {
        if(rank == 0) {
            printf("\n------------------------------------------------\n");
            printf("LAP %d\n\n", lap);
        }
        cur_iter = 0;
        cur_summ = 0;
        null_arr(numbers, MAX_ITERS);
        while(cur_summ < TARGET_SUM) {
            cur_number = generate_number();
            numbers[cur_iter++] = cur_number;
            
            cur_summ += cur_number;
        }

        MPI_Gather(&cur_iter, 1, MPI_INT, iters, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(numbers, MAX_ITERS, MPI_INT, all_numbers, MAX_ITERS, MPI_INT, 0, MPI_COMM_WORLD);

        if(rank == 0) {
            int max_iter = get_max_number(iters, size);
            for(int i = 0; i < max_iter; ++i) {
                printf("iter%d", i + 1);
                for(int j = 0; j < size; ++j) {
                    if(all_numbers[j * MAX_ITERS + i] != 0){
                        printf("    %d", all_numbers[j * 10 + i]);
                        continue;
                    }
                    printf("          ");   
                }
                printf("\n");
            }
            printf("\nRESULTS: ");
            for(int i = 0; i < size; ++i) {
                total_iters[i] += iters[i];
                printf("%d         ", iters[i]);  
            }
            printf("\n------------------------------------------------\n");
        }
    }
    if(rank == 0) {
        printf("TOTAL:   ");
        for(int i = 0; i < size; ++i) {
            printf("%d         ", total_iters[i]);
        }
        printf("\n");
        int min_total_iter = get_min_number(total_iters, size);
        int max_total_iter = get_max_number(total_iters, size);
        printf("Winner:  ");
        for(int i = 0; i < size; ++i) {
            if(total_iters[i] == min_total_iter)
                printf("%d ", i + 1);
        }
        printf("\n");

        printf("Loser:   ");
        for(int i = 0; i < size; ++i) {
            if(total_iters[i] == max_total_iter)
                printf("%d ", i + 1);
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}
