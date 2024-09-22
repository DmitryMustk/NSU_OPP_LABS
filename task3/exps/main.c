#include <stdio.h>
#include <stdlib.h>



int main(void) {
	#pragma omp for
	for (int i = 0; i < 20; ++i) {
		printf("iter %d on thread %d\n", i, );
	}
}
