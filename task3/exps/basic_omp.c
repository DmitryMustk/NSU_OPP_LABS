#include <stdio.h>

int main(void) {
	#pragma omp parallel
	{
		printf("hello world\n");
	}
	return 0;
}
