#!/bin/bash

mpicc -Wall -Wextra main.c && mpiexec -np 4 ./a.out $1 $2 && rm -rf a.out
