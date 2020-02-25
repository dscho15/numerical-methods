# cmake script
cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(NumericalMethods VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(DIR /home/daniel/Desktop/numerical-methods)
set(HEADER ${DIR}/lib)
set(SOURCE1 ${DIR}/src/lecture_1/main.cpp)
set(SOURCE2 ${DIR}/src/lecture_2/main.cpp)

add_executable(lecture1 ${SOURCE1} ${HEADER})
target_include_directories(lecture1 PRIVATE ${HEADER}) 
add_executable(lecture2 ${SOURCE2} ${HEADER})
target_include_directories(lecture2 PRIVATE ${HEADER}) 
