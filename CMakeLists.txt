# cmake script
cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(NumericalMethods VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include_directories(lib)

set(DIR /home/daniel/Desktop/numerical-methods)
set(HEADER ${DIR}/lib)
set(SOURCE1 ${DIR}/src/lecture_1/main.cpp)
set(SOURCE2 ${DIR}/src/lecture_2/main.cpp)
set(SOURCE3 ${DIR}/src/lecture_3/main.cpp)
set(SOURCE4 ${DIR}/src/lecture_4/main.cpp)
set(SOURCE5 ${DIR}/src/mandatory/main.cpp)
set(SOURCE6 ${DIR}/src/lecture_5/main.cpp)
set(SOURCE7 ${DIR}/src/lecture_6/main.cpp)
set(SOURCE8 ${DIR}/src/lecture_7/main.cpp)
set(SOURCE9 ${DIR}/src/lecture_8/main.cpp)


add_executable(lecture1 ${SOURCE1} ${HEADER})
target_include_directories(lecture1 PRIVATE ${HEADER}) 
add_executable(lecture2 ${SOURCE2} ${HEADER})
target_include_directories(lecture2 PRIVATE ${HEADER}) 
add_executable(lecture3 ${SOURCE3} ${HEADER})
target_include_directories(lecture3 PRIVATE ${HEADER}) 
add_executable(lecture4 ${SOURCE4} ${HEADER})
target_include_directories(lecture4 PRIVATE ${HEADER}) 
add_executable(mandatory ${SOURCE5} ${HEADER})
target_include_directories(mandatory PRIVATE ${HEADER}) 
add_executable(lecture5 ${SOURCE6} ${HEADER})
target_include_directories(lecture5 PRIVATE ${HEADER}) 
add_executable(lecture6 ${SOURCE7} ${HEADER})
target_include_directories(lecture6 PRIVATE ${HEADER}) 
add_executable(lecture7 ${SOURCE8} ${HEADER})
target_include_directories(lecture7 PRIVATE ${HEADER}) 
add_executable(lecture8 ${SOURCE9} ${HEADER})
target_include_directories(lecture8 PRIVATE ${HEADER}) 
