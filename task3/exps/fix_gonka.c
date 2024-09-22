#include <stdio.h>

int main(void) {
	int res = 0;
	#pragma omp parallel 
	{
		#pragma omp critical 
		{
			res++;
		}
	}
	printf("res: %d", res);
	return 0;
}
