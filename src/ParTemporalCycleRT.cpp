#include "CycleEnumeration.h"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <ctime>
#include <chrono>
#include <queue>
#include <deque>
#include <list>
#include <set>
#include <iostream>
#include <fstream>

#include <tbb/task.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/spin_mutex.h>
#include <tbb/atomic.h>

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;
using namespace tbb;

extern int timeWindow;
extern bool useCUnion;

namespace {

tbb::combinable<CycleHist> pt_chist_rt;

/// Thread Control block
tbb::atomic<int> globalId = 0;
tbb::enumerable_thread_specific<int> threadIds;

int getTheadId() {
    bool exists = true;
    auto& myid = threadIds.local(exists);
    if(!exists) myid = globalId++;
    return myid;
}

struct ThreadDataGuard {
public:
    ThreadDataGuard() : refCount(1) {
        #ifdef BLK_FORWARD
        blocked = new HashMapStack(true /*concurrent*/);
        #endif
        cycle = new Cycle(true /*concurrent*/);
        tg = new TimestampGroups(true /*concurrent*/);
    }

    ThreadDataGuard(ThreadDataGuard* guard, int lvl, int pathSize) : refCount(1) {
        #ifdef BLK_FORWARD
        blocked = guard->blocked->clone(lvl);
        #endif
        cycle = guard->cycle->clone(pathSize);
        tg = guard->tg->clone(pathSize-1);
    }

    void incrementRefCount() {
        RefCountMutexType::scoped_lock lock(RefCountMutex);
        refCount++;
    }

    void decrementRefCount() {
        RefCountMutexType::scoped_lock lock(RefCountMutex);
        refCount--;
        if(refCount <= 0) {
            #ifdef BLK_FORWARD
            delete blocked; blocked = NULL;
            #endif
            delete cycle; delete tg;
            cycle = NULL; tg = NULL;

            delete this;
        }
    }

    #ifdef BLK_FORWARD
    HashMapStack* blocked = NULL;
    #endif
    Cycle* cycle = NULL;
    TimestampGroups* tg = NULL;

    int refCount = 1;

    typedef tbb::spin_mutex RefCountMutexType;
    RefCountMutexType RefCountMutex;
};

/// Tasks

class TempRTCycleTask : public tbb::task {
public:
    TempRTCycleTask(Graph* gg, int v, EdgeData s, int prevCl, HashMapStack* blk, Cycle* cyc, TimestampGroups* tsg,
            ThreadDataGuard *tdg, TempPathBundle* main, TempPathBundle* alt, StrongComponent* cu = NULL, int lvl = 0):
            g(gg), vert(v), start(s), rootTask(false), thrData(tdg), pathSize(cyc->size()), blocked(blk), cycle(cyc), tg(tsg),
            prevClose(prevCl), level(lvl), main_path(main), alt_path(alt), ind(g->offsArray[v]), cunion(cu)
    { ownerThread = getTheadId();}

    TempRTCycleTask(Graph* gg, int v, int rank = 0, int ncluster = 1):
            g(gg), vert(v), start(v, -1), rootTask(true), level(0), pathSize(0), ind(g->offsArray[v]), process_rank(rank), number_of_clusters(ncluster)
    {
        ind = g->offsArray[v] + (rank + v) % number_of_clusters;
        ownerThread = getTheadId();
    }

    virtual ~TempRTCycleTask(){}
    virtual tbb::task* execute() override;

protected:
    bool isContinuation = false;

    // sequential execution
    tbb::task* outerLoop();
    tbb::task* outerLoopRoot();
    tbb::task* Continuation();

    ///Parameters
    Graph *g = NULL;
    int vert = -1;
    int pathSize;
    EdgeData start = EdgeData(-1, -1);
    int prevClose = -1;

    // Paths
    TempPathBundle* main_path = NULL;
    TempPathBundle* alt_path = NULL;

    ThreadDataGuard *thrData;

    HashMapStack* blocked = NULL;
    Cycle* cycle = NULL;
    TimestampGroups* tg = NULL;

    /// Cycle union
    StrongComponent* cunion = NULL;

    /// MPI control
    int process_rank = 0;
    int number_of_clusters = 1;

    /// Task control
    int seqCnt = 0;
    int level = 0;
    bool rootTask = false;
    int ownerThread = -1;

    // Continuation stealing
    bool newTask = true;
    int ind = 0;

    #ifdef NEIGHBOR_ORDER
    vector<TempPathBundle*> allPaths;
    #endif
    int it = -1;
};

class LoopIterationRT : public tbb::task {
public:
    LoopIterationRT(Graph* gg, int v, TempEdge p, EdgeData s, HashMapStack* blk, Cycle* cyc, TimestampGroups* tsg,
            ThreadDataGuard *tdg, TempPathBundle* cur_path, StrongComponent* cu = NULL, int lvl = 0):
            g(gg), vert(v), prevEdge(p), start(s), blocked(blk), cycle(cyc), tg(tsg), level(lvl), thrData(tdg), cunion(cu),
            current_path(cur_path), rootTask(false) {}

    // Constructor for the root call
    LoopIterationRT(Graph* gg, int v, TempEdge p, EdgeData s):
            g(gg), vert(v), prevEdge(p), start(s), level(0), rootTask(true) {}

    virtual ~LoopIterationRT(){}
    virtual tbb::task* execute() override;

protected:
    bool isContinuation = false;

    tbb::task* innerIteration();
    tbb::task* Continuation();

    ///Parameters
    Graph *g;
    int vert = -1;
    TempEdge prevEdge;
    EdgeData start;

    ThreadDataGuard *thrData;

    HashMapStack* blocked = NULL;
    Cycle* cycle = NULL;
    TimestampGroups* tg = NULL;

    /// Cycle union
    StrongComponent* cunion = NULL;

    /// Paths
    TempPathBundle* current_path = NULL;

    /// Task controll
    int level = 0;
    bool rootTask = false;
};

tbb::task* TempRTCycleTask::execute() {

    if(newTask) {
        set_ref_count(1);
        newTask = false;
    }
    if (!isContinuation) {

        if(!rootTask) {
            // COPY ON STEAL
            int thisThreadId = getTheadId();
            if(ownerThread != thisThreadId) {
                ownerThread = thisThreadId;

                ThreadDataGuard* newb = new ThreadDataGuard(thrData, level, pathSize);
                thrData->decrementRefCount();

                thrData = newb;
                #ifdef BLK_FORWARD
                blocked = thrData->blocked;
                #endif
                cycle = thrData->cycle;
                tg = thrData->tg;
            }
            else
            {
                #ifdef BLK_FORWARD
                blocked->setLevel(level);
                #endif
                cycle->pop_back_until(pathSize);
                tg->pop_back_until(pathSize-1);
            }
        }

        tbb::task* nextTask = nullptr;
        if(rootTask) nextTask = outerLoopRoot();
        else nextTask = outerLoop();

        return nextTask;
    }
    else
        return Continuation();
}

tbb::task* LoopIterationRT::execute() {
    if (!isContinuation) {
        set_ref_count(1);

        tbb::task* nextTask = nullptr;
        nextTask = innerIteration();

        return nextTask;
    }
    else
        return Continuation();
}

tbb::task* TempRTCycleTask::Continuation() {
    return NULL;
}

tbb::task* LoopIterationRT::Continuation() {
    if(rootTask && cunion) {
        delete cunion; cunion = NULL;
    }
    return NULL;
}

/// The main backtracking function
tbb::task* TempRTCycleTask::outerLoop() {
    tbb::task* retTask = NULL;

    /// Iterate through every neighboring vertex and find paths

#if defined(NEIGHBOR_ORDER)
    int ind = g->offsArray[vert];
    if(it == -1)
#endif
    while(ind < g->offsArray[vert+1]) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;
        ind++;

        #if defined(ONE_DAY_CYCLE) || defined(NON_TEMPORAL)
        if (w < startVert) continue;
        #endif

        /// Find the next timestamp that will be used

        // Move index to the starting timestamp where tw >= tg->back().front()
        int prev_min_ts = tg->back().first_ts;
        int end_tstamp = prevClose == 0 ? start.tstamp + timeWindow : prevClose;

        auto it_start = upper_bound(tset.begin(), tset.end(), prev_min_ts);
        auto it_end = upper_bound(it_start, tset.end(), end_tstamp);

        if(it_start >= it_end) continue;

        int tw = *(it_start);

        /// Don't visit this vertex at all
        if ((cunion && !cunion->exists(w)) || blocked->exists(w, tw)) continue;

        /// Find path that needs to be followed
        TempPathBundle *current_path = NULL;
        if (main_path && w == main_path->back().vertex) {
            current_path = main_path;
            main_path = NULL;
        } else if (alt_path && w == alt_path->back().vertex) {
            current_path = alt_path;
            alt_path = NULL;
        } else {
            bool found = false;

            found = findPathTemp(g, TempEdge(w, it_start, it_end), start, *blocked, current_path, cunion);

            if (!found) continue;
        }

#ifdef NEIGHBOR_ORDER
        allPaths.push_back(current_path);
    }

    if(it == -1) it = 0;
    while(it < allPaths.size()) {
        TempPathBundle* current_path = allPaths[it]; it++;
#endif

        // Forwarding the current path
        Cycle* new_cycle = cycle;
        TimestampGroups* new_tg = tg;

        // Forwarding the blocked set
        #ifdef BLK_FORWARD
        HashMapStack* new_blocked = blocked;
        new_blocked->incrementLevel();
        #else
        HashMapStack *new_blocked = new HashMapStack;
        for (auto c : *cycle) if (c != start.vertex) new_blocked->insert(c, -1);
        #endif

        // Creating new task
        increment_ref_count();
        tbb::task *a = new(allocate_child()) LoopIterationRT(g, vert, current_path->back(), start, new_blocked,
                                                             new_cycle, new_tg, thrData, current_path, cunion, level);
#if !defined(NEIGHBOR_ORDER)
        bool exitCond = (ind == g->offsArray[vert + 1]);
#elif defined(NEIGHBOR_ORDER)
        bool exitCond = (it == allPaths.size());
#else
        bool exitCond = (it == allStartingEdges.size());
#endif
        if (!exitCond) {
            recycle_to_reexecute();
            thrData->incrementRefCount();
            return a;
        } else retTask = a;
    }

    if(!retTask) thrData->decrementRefCount();
    #ifndef BLK_FORWARD
    delete blocked; blocked = NULL;
    #endif

    recycle_as_safe_continuation();
    isContinuation = true;

    return retTask;
}

/// The main backtracking function
tbb::task* TempRTCycleTask::outerLoopRoot() {
    tbb::task *retTask = NULL;

    /// Iterate through every neighboring vertex and find paths
    for(int ind = g->offsArray[vert]; ind < g->offsArray[vert+1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;

        for(int j = 0; j < tset.size(); j++) {
            if((ind+j) % number_of_clusters != process_rank) continue;

            TimestampSet::iterator ts_begin = tset.begin() + j;
            TimestampSet::iterator ts_end = tset.begin() + j + 1;

            increment_ref_count();
            tbb::task *a = new(allocate_child()) LoopIterationRT(g, vert, TempEdge(w, ts_begin, ts_end), EdgeData(vert, tset[j]));
            spawn(*a);
        }
    }

    recycle_as_safe_continuation();
    isContinuation = true;

    return retTask;
}

/// The main backtracking function
tbb::task *LoopIterationRT::innerIteration() {
    if (rootTask) {

        //Filter-out vertices not belonging to the same SCC
        if(useCUnion) {
            findCycleUnions(g, EdgeData(prevEdge.vertex, start.tstamp), start.vertex, timeWindow, cunion);
            if(cunion->size() <= 1) {
                delete cunion; cunion = NULL;
                return NULL;
            }
        }

        // Allocate the data structures
        thrData = new ThreadDataGuard();

        #if defined(BLK_FORWARD)
        blocked = thrData->blocked;
        #else
        blocked = new HashMapStack;
        #endif

        bool found = false;
        found = findPathTemp(g, prevEdge, start, *blocked, current_path, cunion);

        if (!found) {
            thrData->decrementRefCount();
            #if !defined(BLK_FORWARD)
            delete blocked;
            #endif
            if(cunion) { delete cunion; cunion = NULL; }
            return NULL;
        }

        cycle = thrData->cycle;
        tg = thrData->tg;

        #ifdef BLK_FORWARD
        blocked->incrementLevel();
        #endif
        cycle->push_back(vert);
    }

    TempPathBundle *another_path = NULL;
    bool branching = false;

    int prev_close = -1;

    /// Traverse the path until we find a branching point
    TempEdge prev;
    while (current_path->back().vertex != start.vertex && !branching) {
        prev = current_path->back();
        current_path->pop_back();

        /// Update the current cycle and path set
        cycle->push_back(prev.vertex);
        prev_close = blocked->at(prev.vertex);
        blocked->insert(prev.vertex, -1);

        int prev_min_ts = *(prev.ts_begin);
        int end_tstamp = prev_close == 0 ? start.tstamp + timeWindow : prev_close;

        /// Extract timestamps and add them to tg
        tg->push_back(TstampInt(prev_min_ts, prev.ts_begin, prev.ts_end));

        /// Detect branching and find an alternate path
        for(int ind = g->offsArray[prev.vertex]; ind < g->offsArray[prev.vertex+1] && !branching; ind++) {
            int u = g->edgeArray[ind].vertex;
            auto& tset = g->edgeArray[ind].tstamps;

            #if defined(ONE_DAY_CYCLE) || defined(NON_TEMPORAL)
            if (u < start.vertex) continue;
            #endif

            // Not an alternative path
            if(u == current_path->back().vertex) continue;

            #if !defined(ONE_DAY_CYCLE) && !defined(NON_TEMPORAL)
            auto it_start = upper_bound(tset.begin(), tset.end(), prev_min_ts);
            #else
            auto it_start = lower_bound(tset.begin(), tset.end(), prev_min_ts);
            #endif

            /// If the timestamp for the given time window exists
            if(it_start == tset.end() || *it_start > end_tstamp) continue;

            int tu = *it_start;

            // If this is the first timestamp after previous and if is a valid alternate path
            if ((!cunion || cunion->exists(u)) && !blocked->exists(u, tu) && !branching) {
                auto it_end = upper_bound(it_start, tset.end(), end_tstamp);

                #ifndef PATH_FORWARD
                branching = dfsPruneTemp(g, TempEdge(u, it_start, it_end), start, *blocked, another_path, cunion, true);
                #else
                branching = dfsPruneTemp(g, TempEdge(u, it_start, it_end), start, *blocked, another_path, cunion, false);
                #endif

                if (branching) break;
            }

        }
    }


    tbb::task *a = NULL;
    /// Record the cycles if needed
    if (!branching) {
        TempEdge finalEdge = current_path->back();
        tg->push_back(TstampInt(*(finalEdge.ts_begin), finalEdge.ts_begin, finalEdge.ts_end));
        recordBundledCycle(cycle, tg, pt_chist_rt.local());
        tg->pop_back();

        #if !defined(BLK_FORWARD)
        delete blocked;
        blocked = NULL;
        #endif

        delete current_path;
        current_path = NULL;

        thrData->decrementRefCount();
    }
    else {
        #ifndef PATH_FORWARD
        if (current_path) delete current_path;
        current_path = NULL;
        #endif

        increment_ref_count();
        a = new(allocate_child()) TempRTCycleTask(g, prev.vertex, start, prev_close, blocked, cycle, tg, thrData, current_path, another_path, cunion, level + 1);
    }

    recycle_as_safe_continuation();
    isContinuation = true;

    return a;
}

/// ******************************* Top level function for root calls *******************************

class RootTempRTCycleTask : public tbb::task {
public:
    RootTempRTCycleTask(Graph *g, int nthr, int rank = 0, int nclust = 1) :
            isContinuation(false),
            graph(g),
            numThreads(nthr),
            process_rank(rank),
            size_of_cluster(nclust) {}

    virtual ~RootTempRTCycleTask() {}

    virtual tbb::task *execute() override;

protected:
    bool isContinuation;
    int numThreads;
    Graph *graph;

    int process_rank;
    int size_of_cluster;
};

tbb::task *RootTempRTCycleTask::execute() {
    if (!isContinuation) {
        set_ref_count(1);

        parallel_for(size_t(0), size_t(graph->getVertexNo()), [&](size_t i) {
            if ((graph->numNeighbors(i) != 0) && (graph->numInEdges(i) != 0)) {
                increment_ref_count();
                spawn(*new(allocate_child()) TempRTCycleTask(graph, i, process_rank, size_of_cluster));
            }
        });

        recycle_as_safe_continuation();
        isContinuation = true;
    }

    return NULL;
}

}

/// ******************************* Fine-grained parallel Read-Tarjan - top level *******************************

void allCyclesTempReadTarjanFineGrained(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    tbb::task *rt = new(RootTempRTCycleTask::allocate_root()) RootTempRTCycleTask(g, numThreads, process_rank, size_of_cluster);
    tbb::task::spawn_root_and_wait(*rt);

    combineCycleHistogram(pt_chist_rt, result);
}