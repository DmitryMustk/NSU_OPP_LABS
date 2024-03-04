cmake_minimum_required(VERSION 3.10)

project(task1)

# Определение переменных SOURCES
set(SOURCES)

# Версия 1
if(VERSION STREQUAL "v1")
    message("Building v1 (sequential programm)")
    list(APPEND SOURCES
        src/main.c
    )
    add_executable(task1 ${SOURCES})
    set(CMAKE_C_COMPILER "gcc")
endif()

# Версия 2
if(VERSION STREQUAL "v2")
    message("Building v2 (paralell with dot to dot communication)")
    list(APPEND SOURCES
        src/main2.c
    )
    add_executable(task1 ${SOURCES})

    set(CMAKE_C_COMPILER "mpicc")
endif()

# Версия 3
if(VERSION STREQUAL "v3")
    message("Building v3 (parallel with collective communication)")
    list(APPEND SOURCES
        src/main3.c
    )
    add_executable(task1 ${SOURCES})
    set(CMAKE_C_COMPILER "mpicc")
endif()
