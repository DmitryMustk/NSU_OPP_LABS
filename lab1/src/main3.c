#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 100000

int main(int argc, char **argv) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start, end;
    int vec_part_size = N / size;

    double* a = malloc(sizeof(double) * N);
    double* b = rank == 0 ? malloc(sizeof(double) * N) : NULL; 
    double* rec = malloc(sizeof(double) * vec_part_size);           

    if (rank == 0){
        for (int i = 0; i < N; ++i){
            a[i] = 1.0;
            b[i] = 1.0;
        }

        start = MPI_Wtime();
    }

    MPI_Bcast(a, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, vec_part_size, MPI_DOUBLE, rec, vec_part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double part_sum = 0;
    for (int j = 0; j < vec_part_size; ++j){
        for (int i = 0; i < N; ++i)
            part_sum += a[i] * rec[j];
    }

    if (rank == 0){
        for (int j = vec_part_size * size; j < N; ++j){
            for (int i = 0; i < N; ++i)
                part_sum += a[i] * b[j];
        }
    }

    double s = 0;
    MPI_Reduce(&part_sum, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        end = MPI_Wtime();
        printf("Time taken: %f\n", end - start);
        printf("Result: %f\n", s);
    }

    free(a);
    free(b);
    free(rec);

    MPI_Finalize();
    return 0;
}
