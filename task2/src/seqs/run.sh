#!/bin/bash

set -xe

mpicc -Wall -Wextra $1 -lm && mpiexec -np $2 a.out && rm -rf a.out
