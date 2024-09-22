#!/bin/bash

mpicxx ../exps/main.cpp && mpiexec -np $2 ./a.out $3 $4 && rm -rf a.out
