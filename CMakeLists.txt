cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

project(vices VERSION 1.0.0)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(ZLIB REQUIRED)
find_library(STATGEN_LIBRARY StatGen)
find_package(Armadillo CONFIG REQUIRED)

add_executable(vices src/main.cpp src/regression.cpp src/readReports.cpp src/SampleData.cpp src/Report.cpp src/VicesArgs.cpp)

target_link_libraries(vices ${STATGEN_LIBRARY} ${ZLIB_LIBRARY} armadillo)

install(TARGETS vices DESTINATION bin)
