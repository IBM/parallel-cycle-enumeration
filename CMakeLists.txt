cmake_minimum_required(VERSION 3.2)
project(parallel_cycle_enumeration)

INCLUDE($ENV{TBBROOT}/cmake/TBBBuild.cmake)
SET(TBB_DIR "$ENV{TBBROOT}/cmake")

SET(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_COLOR_MAKEFILE ON)

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add here all include directories
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

SET(SRC_FILES
        src/CycleEnumeration.cpp
        src/CycleReadTarjan.cpp
        src/JohnsonsAlgorithm.cpp
        src/TempCycleJohnson.cpp
        src/TemporalCycleRT.cpp
        src/ParTemporalCycleRT.cpp
        src/ParCycleRT.cpp
        src/CycleMain.cpp)

ADD_EXECUTABLE(cycle ${SRC_FILES})

FIND_PACKAGE(TBB COMPONENTS tbbmalloc tbbmalloc_proxy REQUIRED)
TARGET_LINK_LIBRARIES(cycle tbb tbbmalloc tbbmalloc_proxy)
INCLUDE_DIRECTORIES(${TBB_DIR}/../include)