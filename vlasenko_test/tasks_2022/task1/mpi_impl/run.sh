#!/bin/bash

./build.sh && cd build && mpirun -np 6 ./2022_task1_mpi_impl 6 
