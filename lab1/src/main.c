#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <time.h>

#define N 100000

void fill_vec(double* vec){
    for(size_t i = 0; i < N; ++i){
        vec[i] = 1.0;
    }
}

double calculate(double* a, double* b){
    double s = 0;
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            s += a[i] * b[j];
        }
    }
}

int main(void) {
    srand(time(NULL));
    double* a = malloc(sizeof(double) * N);
    double* b = malloc(sizeof(double) * N);

    fill_vec(a);
    fill_vec(b);
    clock_t start = clock();
    int64_t res = calculate(a, b);
    clock_t end = clock();
    printf("Time taken: %.2f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Result: %ld\n", res);
    free(a);
    free(b);
    return 0;
}
