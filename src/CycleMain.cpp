#include <iostream>
#include <fstream>

#include <cstdint>
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>
#include <vector>
#include <ctime>
#include <chrono>
#include <map>

#include "utils.h"
#include "Graph.h"
#include "CycleEnumeration.h"

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;

int timeWindow = 3600;
bool temporal = true;
bool useCUnion = false;
ConcurrentCounter vertexVisits;

// *************************** function main ***********************************

int main(int argc, char** argv) {
    int process_rank = 0; int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    if(argc == 1 || cmdOptionExists(argv, argv+argc, "-h")) {
        printHelp();
        return 0;
    }

    string path;
    if(cmdOptionExists(argv, argv+argc, "-f")) path = string(getCmdOption(argv, argv + argc, "-f"));
    else cout << "No input file" << endl;
    string name = path.substr(path.find_last_of("/") + 1, path.find_last_of(".")-path.find_last_of("/")-1);

    /// Check if the file exists
    if(!pathExists(path)) { cout << "The input file doesn't exist" << endl; return 0; }

    int nthr = 256;
    int algo = -1;
    if(cmdOptionExists(argv, argv+argc, "-n")) nthr = stol(string(getCmdOption(argv, argv + argc, "-n")));
    if(cmdOptionExists(argv, argv+argc, "-tw")) timeWindow = stol(string(getCmdOption(argv, argv + argc, "-tw")));
    if(cmdOptionExists(argv, argv+argc, "-tws")) timeWindow = stol(string(getCmdOption(argv, argv + argc, "-tws")));
    if(cmdOptionExists(argv, argv+argc, "-algo")) algo = stoi(string(getCmdOption(argv, argv + argc, "-algo")));
    if(cmdOptionExists(argv, argv+argc, "-cunion")) useCUnion = true;
    if(cmdOptionExists(argv, argv+argc, "-tw")) timeWindow *= 3600;

    if(process_rank == 0) cout << "Reading " << name << endl;
    auto read_start = chrono::steady_clock::now();
    Graph *g = new Graph; g->readTemporalGraph(path);
    auto read_end = chrono::steady_clock::now();
    double read_total = chrono::duration_cast<chrono::milliseconds>(read_end - read_start).count()/1000.0;
    if(process_rank == 0) cout << "Graph read time time: " << read_total << " s" << endl;

    #ifdef MPI_IMPL
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    bool printResult = true;

    CycleHist resultHistogram;
    auto total_start = chrono::steady_clock::now();
    switch(algo) {

        case 0:
            if(process_rank == 0) cout << " -----------------  Coarse-grained parallel JOHNSON with time window ------------------ " << endl;
            allCyclesJohnsonCoarseGrainedTW(g, resultHistogram, nthr);
            break;

        case 1:
            if(process_rank == 0) cout << " ----------------- Fine-grained parallel JOHNSON with time window ------------------ " << endl;
            allCyclesJohnsonFineGrainedTW(g, resultHistogram, nthr);
            break;

        case 2:
            if(process_rank == 0) cout << " ----------------- Coarse-grained parallel READ-TARJAN with time window ------------------ " << endl;
            allCyclesReadTarjanCoarseGrainedTW(g, resultHistogram, nthr);
            break;

        case 3:
            if(process_rank == 0) cout << " ----------------- Fine-grained parallel READ-TARJAN  with time window ------------------ " << endl;
            allCyclesReadTarjanFineGrainedTW(g, resultHistogram, nthr);
            break;

        case 4:
            if(process_rank == 0)  cout << " ----------------- Coarse-grained parallel temporal JOHNSON ----------------- " << endl;
            allCyclesTempJohnsonCoarseGrained(g, resultHistogram, nthr); break;

        case 5:
            if(process_rank == 0) cout << " ----------------- Fine-grained parallel temporal JOHNSON ------------------ " << endl;
            allCyclesTempJohnsonFineGrained(g, resultHistogram, nthr);
            break;

        case 6:
            if(process_rank == 0) cout << " ----------------- Coarse-grained parallel temporal READ-TARJAN ------------------ " << endl;
            allCyclesTempReadTarjanCoarseGrained(g, resultHistogram, nthr);
            break;

        case 7:
            if(process_rank == 0) cout << " ----------------- Fine-grained parallel temporal READ-TARJAN ------------------ " << endl;
            allCyclesTempReadTarjanFineGrained(g, resultHistogram, nthr);
            break;

        case 8:
            cout << " ----------------- Temporal Read-Tarjan ------------------ " << endl;
            allCyclesReadTarjanTemp(g, resultHistogram);
            break;

        case 9:
            cout << " ------------ 2SCENT without source detection  ------------ " << endl;
            allCycles2scent(g, resultHistogram);
            break;

        case 10:
            cout << " ------------ 2SCENT with source detection ------------ " << endl;
            {
                SeedSet* seedset = NULL;
                double seed_tot = createSeeds2scent(g, seedset);
                cout << "Preprocessing time: " << seed_tot << " s" << endl;

                auto total_start = chrono::steady_clock::now();
                allCycles2scentSeed(g, *seedset, resultHistogram);
                auto total_end = chrono::steady_clock::now();
                double total = chrono::duration_cast<chrono::milliseconds>(total_end - total_start).count()/1000.0;
                cout << "Backtracking time: " << total << " s" << endl;

                delete seedset;
            }
            break;

        case 11:
            cout << " ------------ Cycle-union execution time ------------ " << endl;
            cycleUnionExecTime(g, nthr);
            break;

        default:
            printHelp();
            printResult = false;
            break;
    }

#ifdef MPI_IMPL
    if (process_rank == 0) {
            for(int r = 1; r < size_of_cluster; r++) {
                unsigned long stat;
                MPI_Recv(&stat, 1, MPI_UNSIGNED_LONG, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                vertexVisits.increment(stat);
            }
        } else {
            unsigned long stat = vertexVisits.getResult();
            MPI_Send(&stat, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        }
#endif

    if(process_rank == 0) {
        auto total_end = chrono::steady_clock::now();
        double total = chrono::duration_cast<chrono::milliseconds>(total_end - total_start).count()/1000.0;
        cout << "Total time: " << total << " s" << endl;
        cout << "Vertex visits: " << vertexVisits.getResult() << endl;
        if (printResult) {
            cout << "# cycle_size, num_of_elems\n";
            unsigned long totCycles = 0;
            for (auto hist : resultHistogram) { cout << hist.first << ", " << hist.second << "\n"; totCycles += hist.second; }
            cout << "Total, " << totCycles << endl;
        }
    }
    delete  g;

#ifdef MPI_IMPL
    MPI_Finalize();
#endif

    return 0;
}
