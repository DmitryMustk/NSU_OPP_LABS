#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 100000

int main(int argc, char **argv){
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vec_part_size = N / size; 

    if (rank == 0){
        double* a = malloc(sizeof(double) * N);
        double* b = malloc(sizeof(double) * N);
        for (int i = 0; i < N; ++i){
            a[i] = 1.0;
            b[i] = 1.0;
        }
        double s = 0;
        double start = MPI_Wtime();

        for (int p = 1; p < size; ++p){
            MPI_Send(a, N, MPI_DOUBLE, p, 1, MPI_COMM_WORLD);
            MPI_Send(b + vec_part_size * (p - 1), vec_part_size, MPI_DOUBLE, p, 2, MPI_COMM_WORLD);
        }

        for (int j = vec_part_size * (size - 1); j < N; ++j){
            for (int i = 0; i < N; ++i)
                s += a[i] * b[j];
        }

        MPI_Status st;
        double part_sum;
        for (int p = 1; p < size; ++p) {
            //MPI_Probe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &st);
            MPI_Recv(&part_sum, 1, MPI_DOUBLE, st.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            s += part_sum;
        }

        double end = MPI_Wtime();
        printf("Time taken: %f\n", end - start);
        printf("Result: %f\n", s);

        free(a);
        free(b);
    }
    else{
        double* a = malloc(sizeof(double) * N);
        double* b = malloc(sizeof(double) * vec_part_size);

        MPI_Recv(a, N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, vec_part_size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double part_sum = 0;
        for (int j = 0; j < vec_part_size; ++j){
            for (int i = 0; i < N; ++i)
                part_sum += a[i] * b[j];
        }

        MPI_Send(&part_sum, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

        free(a);
        free(b);
    }

    MPI_Finalize();
    return 0;
}