cmake_minimum_required(VERSION 3.0)
project(2022_task1_mpi_impl LANGUAGES C)

find_package(MPI REQUIRED)

add_executable(${PROJECT_NAME} src/main.c)
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_C)

