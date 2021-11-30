# parallel_cycle_enumeration

Source code for parallel enumeration of simple and temporal cycles.

The contact email address: jov@zurich.ibm.com

### Building

Prerequisites for building our code:
- GCC version 7 or higher
- MPICC compiler
- Intel Threading Building Blocks 2020 Update 2, available here: https://github.com/oneapi-src/oneTBB/releases/tag/v2020.2

Before building our code, tbb needs to be enabled:

```
source /path-to-tbb/bin/tbbvars.sh intel64 linux auto_tbbroot
```

If MPI is available, built this code using the following commands:

```
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_CXX_COMPILER=`which mpicxx`
make
```

To run it on multiple nodes, comment out `MPI_IMPL` macro from `include/Macros.h` file and use:
```
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`
make
```

### Executing

To run our code on a single node, simply execute:
```
./cycle -f <graph_path> -algo <algo> -tw <time-window> -n <num-of-threads>
```
Or using `mpirun`:

```
mpirun -np <num-of-processes> --ppn 1 -hosts <list-of-hosts> ./cycle -f <graph_path> -algo <algo> -tw <time-window> -n <num-of-threads>
```

Optionally, for better performance, interleave the memory across different NUMA nodes using numactl command:
```
numactl -i all ./cycle <cmd-arguments>
```

When running the code using `mpicc` on multiple processors, make sure to have a copy of the input graph in the same path across all nodes.
Also, the path to `tbb` and the path to the executable should be the same across all nodes.

The exact algorithms supported and the description of the command line arguments can be printed using `./cycle -h`: 
```
 Parallel cycle enumeration algorithms 
    -f             Path to the input file graph
    -n             Number of threads, default 256
    -tw            Time window in hours, default 1h
    -tws           Time window in seconds, default 1h
    -cunion        Use cycle-unions for temporal cycle enumeration
    -algo          Define which algorithm to execute
                       0  - Coarse-grained parallel JOHNSON with time window
                       1  - Fine-grained parallel JOHNSON with time window
                       2  - Coarse-grained parallel READ-TARJAN with time window
                       3  - Fine-grained parallel READ-TARJAN  with time window
                       4  - Coarse-grained parallel temporal JOHNSON
                       5  - Fine-grained parallel temporal JOHNSON
                       6  - Coarse-grained parallel temporal READ-TARJAN
                       7  - Fine-grained parallel temporal READ-TARJAN
                       8  - Temporal Read-Tarjan - single thread
                       9  - 2SCENT without source detection - single thread
                       10 - 2SCENT with source detection - single thread
                       11 - Cycle-union execution time test
    -h             Print this message
```

The file containing the input temporal graph contains the list of edges, one edge per line. Each edge is represented with the identifiers of the two vertices that it connects and a timestamp.
Values in each line are separated by space.

The sample input graph given in `data/sample.txt`:

```
1 2 0
1 5 1
2 3 2
3 1 3
3 2 4
3 4 5
3 6 6
4 5 7
5 2 8
6 4 9
4 7 10
```

An example of running our code:

```
numactl -i all ./cycle -f ../data/sample.txt -algo 0 -n 16 -tw 1
```

which outputs the following cycle histogram:
```
# cycle_size, num_of_cycles
2, 1
3, 1
4, 2
5, 1
Total, 5
```