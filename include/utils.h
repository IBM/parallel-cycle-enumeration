#ifndef _GRAPH_UTILS_H_
#define _GRAPH_UTILS_H_

#include <cstdint>
#include <sstream>
#include <fstream>
#include <tbb/combinable.h>
#include <tbb/spin_mutex.h>
#include <tbb/atomic.h>
#include <sys/mman.h>
#include <iostream>
#include <sys/stat.h>
#include <algorithm>
#include <string>

#include "Macros.h"

using namespace std;

inline bool pathExists(string path) {
    struct stat buffer;
    if (stat(path.c_str(), &buffer) != 0) return false;
    return true;
}

inline char* getCmdOption(char ** begin, char ** end, const std::string & option) {
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) return *itr;
    return 0;
}

inline bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

inline void printHelp()
{
    std::cout << " Parallel cycle enumeration algorithms " << std::endl;
    std::cout << "    -f             Path to the input file graph" << std::endl;
    std::cout << "    -n             Number of threads, default 256" << std::endl;
    std::cout << "    -tw            Time window in hours, default 1h" << std::endl;
    std::cout << "    -tws           Time window in seconds, default 1h" << std::endl;
    std::cout << "    -cunion        Use cycle-unions for temporal cycle enumeration" << std::endl;
    std::cout << "    -algo          Define which algorithm to execute" << std::endl;
    std::cout << "                       0  - Coarse-grained parallel JOHNSON with time window" << std::endl;
    std::cout << "                       1  - Fine-grained parallel JOHNSON with time window" << std::endl;
    std::cout << "                       2  - Coarse-grained parallel READ-TARJAN with time window" << std::endl;
    std::cout << "                       3  - Fine-grained parallel READ-TARJAN  with time window" << std::endl;
    std::cout << "                       4  - Coarse-grained parallel temporal JOHNSON" << std::endl;
    std::cout << "                       5  - Fine-grained parallel temporal JOHNSON" << std::endl;
    std::cout << "                       6  - Coarse-grained parallel temporal READ-TARJAN" << std::endl;
    std::cout << "                       7  - Fine-grained parallel temporal READ-TARJAN" << std::endl;
    std::cout << "                       8  - Temporal Read-Tarjan - single thread" << std::endl;
    std::cout << "                       9  - 2SCENT without source detection - single thread" << std::endl;
    std::cout << "                       10 - 2SCENT with source detection - single thread" << std::endl;
    std::cout << "                       11 - Cycle-union execution time test" << std::endl;
    std::cout << "    -h             Print this message" << std::endl;
}

#endif//_GRAPH_UTILS_H_