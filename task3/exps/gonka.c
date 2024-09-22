#include <stdio.h>

int main(void) {
	int res = 0;
	#pragma omp parallel
	{
		res++;
	}
	printf("res: %d\n", res);
	return 0;
}
