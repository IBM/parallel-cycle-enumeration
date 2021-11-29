#include "CycleEnumeration.h"
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <utility>
#include <ctime>
#include <chrono>
#include <queue>
#include <tbb/spin_mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;
using namespace tbb;

extern int timeWindow;

bool edgeInTimeInterval(int tstart, int timeWindow, int vstart, int vert, TimestampSet& tset) {
    auto it_start = (vstart > vert) ? lower_bound(tset.begin(), tset.end(), tstart) :
                    upper_bound(tset.begin(), tset.end(), tstart);
    auto it_end = upper_bound(it_start, tset.end(), tstart + timeWindow);
    if (it_start >= it_end) return false;
    return true;
}

/// *************************** Recording cycles ***************************

void recordCycle(Cycle* current, CycleHist& result) {
    int size = current->size();
    if(result.find(size) == result.end()) result[size] = 0;
    result[size]++;
}

void recordBundledCycle(Cycle* current, TimestampGroups* tg, CycleHist& result) {
    /// Check if the budle is consisted of exactly one cycle
    bool allone = true;
    long count = 1;

    for (auto it = tg->begin(); it != tg->end(); ++it) {
        if (it->it_end - it->it_begin > 1) {
            allone = false;
            break;
        }
    }

    /// Count the number of cycles in the bundle (2SCENT procedure)
    if (!allone) {
        queue <pair<int, long>> *prevQueue = new queue <pair<int, long>>;
        queue <pair<int, long>> *currQueue = nullptr;

        prevQueue->push(make_pair(-1, 1));

        for (auto it = tg->begin(); it != tg->end(); ++it) {
            currQueue = new queue <pair<int, long>>;
            int n = 0, prev = 0;
            for (auto it2 = it->it_begin; it2 != it->it_end; ++it2) {
                int ts = *it2;

                if (!prevQueue->empty()) {
                    auto tmpPair = prevQueue->front();
                    while (tmpPair.first < ts) {
                        prevQueue->pop(); n = tmpPair.second;
                        if (prevQueue->empty()) break;
                        tmpPair = prevQueue->front();
                    }
                }
                prev += n;
                currQueue->push(make_pair(ts, prev));
            }
            delete prevQueue;
            prevQueue = currQueue;
            currQueue = nullptr;
        }
        auto tmpPair = prevQueue->back();
        count = tmpPair.second;
        delete prevQueue;
    }

    long size = current->size();
    if (result.find(size) == result.end()) result[size] = 0;
    result[size] += count;
}

void combineCycleHistogram(tbb::combinable<CycleHist>& pt_chist, CycleHist& result) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    CycleHist histogram;

    pt_chist.combine_each([&](CycleHist hist) {
        for (auto &pair : hist) {
            const int &cyc_size = pair.first;
            const unsigned long &cyc_num = pair.second;

            if (histogram.find(cyc_size) == histogram.end())
                histogram[cyc_size] = 0;

            histogram[cyc_size] += cyc_num;
        }
    });
#ifdef MPI_IMPL
    // Tag: 0 - hist size; 1 - histdata
    if (process_rank == 0) {
        for(int r = 1; r < size_of_cluster; r++) {
            int histsize = -1;
            MPI_Recv(&histsize, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            long* rec_hist = new long [2*histsize];
            MPI_Recv(rec_hist, 2*histsize, MPI_LONG, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int i = 0; i < histsize; i++) {
                if (histogram.find(rec_hist[2*i]) == histogram.end()) histogram[rec_hist[2*i]] = 0;
                histogram[rec_hist[2*i]] += rec_hist[2*i+1];
            }

            delete [] rec_hist;
        }

    } else{
        int histsize = histogram.size();
        MPI_Send(&histsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        long* send_hist = new long [2*histsize];

        int i = 0;
        for (auto hist : histogram) { send_hist[i++] = hist.first; send_hist[i++] = hist.second; }

        MPI_Send(send_hist, 2*histsize, MPI_LONG, 0, 1, MPI_COMM_WORLD);

        delete [] send_hist;
    }
#endif

    if(process_rank == 0) result = histogram;
}

/// *************************** Find Cycle-Unions ***********************************

const int PINF = std::numeric_limits<int>::max();

void findTempAncestors(Graph *g, EdgeData u, int firstTs, int timeWindow, BlockedMap& visited, BlockedMap* candidates = NULL) {
    /// Temporary prevent re-visiting this vertex
    visited.insert(u.vertex, PINF);

    int mints = PINF;
    int end_tstamp = firstTs + timeWindow;

    for(int ind = g->inOffsArray[u.vertex + 1] - 1; ind >= g->inOffsArray[u.vertex]; ind--) {
        int w = g->inEdgeArray[ind].vertex;
        auto& tset = g->inEdgeArray[ind].tstamps;

        auto it_start = lower_bound(tset.begin(), tset.end(), u.tstamp);

        /// Update maxts
        if(it_start != tset.end()) {
            mints = min(mints, *it_start);
        }

        /// If the timestamp is within the time interval
        if(it_start == tset.begin() || *(it_start-1) <= firstTs) continue;

        int ts = *(it_start-1);

        if(candidates && !candidates->exists(w)) continue;

        /// Recursively visit the neighbor if it can be visited
        if (!visited.exists(w) || ts > visited.at(w)) {
            findTempAncestors(g, EdgeData(w, ts), firstTs, timeWindow, visited, candidates);
        }
    }
    /// Update the closing time value
    visited.insert(u.vertex, mints);
}

int findCycleUnions(Graph *g, EdgeData startEdge, int startVert, int timeWindow, StrongComponent*& cunion) {
    BlockedMap tempAncestors;
    tempAncestors.insert(startVert, PINF);

    int firstTs = startEdge.tstamp;
    for(int ind = g->inOffsArray[startVert + 1] - 1; ind >= g->inOffsArray[startVert]; ind--) {
        int w = g->inEdgeArray[ind].vertex;
        auto& tset = g->inEdgeArray[ind].tstamps;

        auto it_end = upper_bound(tset.begin(), tset.end(), firstTs + timeWindow);
        if(it_end == tset.begin() || *(it_end-1) <= firstTs) continue;

        findTempAncestors(g, EdgeData(w, *(it_end-1)), firstTs,  timeWindow, tempAncestors); // , &tempDescendants
    }

    cunion = new StrongComponent(g->getVertexNo());
    tempAncestors.for_each([&](int el){ cunion->insert(el); });

    return 0;
}

void cycleUnionExecTime(Graph *g, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    parallel_for(size_t(0), size_t(g->getVertexNo()), [&](size_t vert) {
        parallel_for(size_t(g->offsArray[vert]), size_t(g->offsArray[vert + 1]), [&](size_t ind) {
            int w = g->edgeArray[ind].vertex;
            auto &tset = g->edgeArray[ind].tstamps;

            parallel_for(size_t(0), size_t(tset.size()), [&](size_t j) {

                if ((ind + j) % size_of_cluster == process_rank) {
                    int tw = tset[j];
                    StrongComponent *cunion = NULL;
                    findCycleUnions(g, EdgeData(w, tw), vert, timeWindow, cunion);
                    delete cunion;
                }
            });
        });
    });
}