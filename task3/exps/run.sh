#/bin/bash

gcc -fopenmp -Wall -Wextra -lm $1 && ./a.out && rm -rf a.out

