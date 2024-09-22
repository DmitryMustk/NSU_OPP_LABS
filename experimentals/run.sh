#!/bin/bash

mpicc $1 -lm && mpiexec -np $2 a.out && rm -rf a.out
