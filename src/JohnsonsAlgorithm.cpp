#include "CycleEnumeration.h"
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <stack>
#include <utility>
#include <ctime>
#include <chrono>
#include <queue>
#include <vector>


#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;
using namespace tbb;

extern int timeWindow;
extern ConcurrentCounter vertexVisits;

namespace {

tbb::combinable <CycleHist> pt_chist_jh;

/// Blocking and unblocking vertices

void unblockJohnson(int node, BlockedSet &blocked, BlockedList &BMap) {
    blocked.remove(node);
    auto it = BMap.find(node);
    if (it == BMap.end()) return;
    it->second.for_each([&](int w) {
        if (blocked.exists(w)) {
            unblockJohnson(w, blocked, BMap);
        }
    });
    BMap.erase(it);
}

void addToBMap(Graph *g, int node, BlockedList &BMap, StrongComponent *scc = NULL, int s = -1, int tstart = -1) {
    for (int ind = g->offsArray[node]; ind < g->offsArray[node + 1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto &tset = g->edgeArray[ind].tstamps;
        if (s != -1 && tstart != -1 && !edgeInTimeInterval(tstart, timeWindow, s, node, tset)) continue;
        if (!scc || scc->exists(w)) {
            BMap[w].insert(node);
        }
    }
}

/// Thread Control block
tbb::atomic<int> globalId = 0;
tbb::enumerable_thread_specific<int> threadIds;

int getThreadId() {
    bool exists = true;
    auto &myid = threadIds.local(exists);
    if (!exists) myid = globalId++;
    return myid;
}

typedef tbb::mutex ThrDataMutexType;
typedef tbb::spin_mutex RefCountMutexType;

/// ********************** Task-parallel Johnson's algorithm ***********************

struct ThreadDataGuard {
public:
    ThreadDataGuard(Cycle *c, BlockedSet *b, BlockedList *bm) :
            cycle(c), blocked(b), blkMap(bm), refCount(1) {}

    ThreadDataGuard(ThreadDataGuard *guard, int pathSize, int ind) : refCount(1) {
        // Copy the data from another thread
        {
            ThrDataMutexType::scoped_lock lock(guard->ThrDataMutex);

            blocked = new BlockedSet(*(guard->blocked));
            blkMap = new BlockedList(*(guard->blkMap));
            cycle = new Cycle(*(guard->cycle));
        }

        // Remove the invalid vertices
        while (cycle->size() > pathSize) {
            int lastVertex = cycle->back();
            cycle->pop_back();
            unblockJohnson(lastVertex, *blocked, *blkMap);
        }
    }

    void incrementRefCount() {
        RefCountMutexType::scoped_lock lock(RefCountMutex);
        refCount++;
    }

    void decrementRefCount() {
        RefCountMutexType::scoped_lock lock(RefCountMutex);
        refCount--;
        if (refCount <= 0) {
            delete cycle;
            delete blocked;
            delete blkMap;
            cycle = NULL;
            blocked = NULL;
            blkMap = NULL;
            delete this;
        }
    }

    /// Guarded data
    Cycle *cycle = NULL;
    BlockedSet *blocked = NULL;
    BlockedList *blkMap = NULL;

    ThrDataMutexType ThrDataMutex;
private:
    // Reference counter
    int refCount = 1;

    RefCountMutexType RefCountMutex;
};

class JohnsonsTask : public tbb::task {
public:
    JohnsonsTask(Graph *g, int v, Cycle *cyc, BlockedSet *blk, BlockedList *bmap,
                 ThreadDataGuard *tdg, JohnsonsTask *par = NULL) :
            graph(g), vert(v), parent(par), cycle(cyc), blocked(blk), pathSize(cyc->size()), blkMap(bmap),
            ownerThread(getThreadId()), thrData(tdg), ind(g->offsArray[vert]) {}

    JohnsonsTask(Graph *g, int v, int tst, Cycle *cyc, BlockedSet *blk, BlockedList *bmap,
                 ThreadDataGuard *tdg, JohnsonsTask *par = NULL) :
            graph(g), vert(v), parent(par), cycle(cyc), tstart(tst), blocked(blk), pathSize(cyc->size()),
            blkMap(bmap), ind(g->offsArray[vert]), ownerThread(getThreadId()), thrData(tdg)
            {}

    virtual ~JohnsonsTask() {}

    virtual tbb::task *execute() override;

    void returnValue(bool found, int thrId);
//protected:

    tbb::task *SpawnTask();

    tbb::task *Continuation();

    //Parameters
    int vert;
    int tstart = -1;
    Cycle *cycle = NULL;
    BlockedSet *blocked = NULL;
    BlockedList *blkMap = NULL;
    Graph *graph = NULL;

    //Return
    JohnsonsTask *parent = NULL;
    tbb::atomic<bool> retCycleFound = false;

    // Continuation stealing
    bool newTask = true;
    int ind = -1;

    // Task control
    bool isContinuation = false;
    int ownerThread = -1;
    bool stolenTask = false;

    int pathSize = 0;
    ThreadDataGuard *thrData = NULL;
};

tbb::task *JohnsonsTask::execute() {

    tbb::task *ret = NULL;
    if (!isContinuation) ret = SpawnTask();
    else ret = Continuation();

    return ret;
}

tbb::task *JohnsonsTask::SpawnTask() {
    /// Copy-on-steal
    int thisThreadId = getThreadId();
    if (ownerThread != thisThreadId) {
        stolenTask = true;
        ownerThread = thisThreadId;

        // Copy on steal
        ThreadDataGuard *newThrData = new ThreadDataGuard(thrData, pathSize, vert);

        // Decrement the ref. count of the previous blocked map
        thrData->decrementRefCount();

        thrData = newThrData;
        // Update the pointers
        cycle = thrData->cycle;
        blocked = thrData->blocked;
        blkMap = thrData->blkMap;
    }

    /// If we are executing the task for the first time
    if (newTask) {
        set_ref_count(1);
        newTask = false;
        vertexVisits.increment();
        {
            ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
            cycle->push_back(vert);
            pathSize++;
            blocked->insert(vert);
        }
    }

    tbb::task *retTask = NULL;

    /// Iterate through the neighbors
    while (ind < graph->offsArray[vert + 1]) {
        int w = graph->edgeArray[ind].vertex;
        auto &tset = graph->edgeArray[ind].tstamps;
        ind++;

        if (tstart != -1) {
            if (!edgeInTimeInterval(tstart, timeWindow, cycle->front(), vert, tset)) continue;
        }

        if (w == cycle->front()) {
            auto &my_hist = pt_chist_jh.local();
            recordCycle(cycle, my_hist);
            retCycleFound = true;
        } else if (((tstart != -1) || (tstart == -1) && (w > cycle->front())) && !blocked->exists(w)) {

           JohnsonsTask *a = new(allocate_child()) JohnsonsTask(graph, w, tstart, cycle, blocked, blkMap, thrData, this);

            increment_ref_count();
            /// Continuation stealing
            thrData->incrementRefCount();
            if (ind < graph->offsArray[vert + 1]) {
                recycle_to_reexecute();
                return a;
            } else retTask = a;
        }
    }

    // No continuation happens if the task is stolen
    if (stolenTask) thrData->decrementRefCount();

    recycle_as_safe_continuation();
    isContinuation = true;

    return retTask;
}

void JohnsonsTask::returnValue(bool found, int thrId) {
    if (found) retCycleFound = true;
};

tbb::task *JohnsonsTask::Continuation() {
    if (!stolenTask) {

        if (getThreadId() == ownerThread && pathSize == cycle->size()) {
            int v = cycle->back();

            if (retCycleFound) {
                ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
                unblockJohnson(v, *thrData->blocked, *thrData->blkMap);
            } else {
                ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
                if (tstart == -1) addToBMap(graph, v, *thrData->blkMap, NULL);
                else addToBMap(graph, v, *thrData->blkMap, NULL, cycle->front(), tstart);
            }

            cycle->pop_back();

            // Return
            if (parent) parent->returnValue(retCycleFound, ownerThread);
        }

        thrData->decrementRefCount();
    }

    return NULL;
}
}
/// ************ Johnson algorithm with time window ************

// Johnsons algorithm
bool cyclesJohnsonTW(Graph* g, int vert, int tstart, Cycle* cycle, BlockedSet& blocked,
                     BlockedList&  BMap, CycleHist& result) {

    bool cycleFound = false;
    cycle->push_back(vert);
    blocked.insert(vert);

    vertexVisits.increment();

    for(int ind = g->offsArray[vert]; ind < g->offsArray[vert+1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;

        if(!edgeInTimeInterval(tstart, timeWindow, cycle->front(), vert, tset)) continue;

        if (w == cycle->front()) {
            recordCycle(cycle, result);
            cycleFound = true;
        } else if (!blocked.exists(w)) {
            if (cyclesJohnsonTW(g, w, tstart, cycle, blocked, BMap, result)) cycleFound = true;
        }
    }

    cycle->pop_back();
    if (cycleFound) unblockJohnson(vert, blocked, BMap);
    else addToBMap(g, vert, BMap, NULL, cycle->front(), tstart);

    return cycleFound;
}

namespace {
class RootJohnTask : public tbb::task {
public:
    RootJohnTask(Graph *_g, int nthr, int rank = 0, int nclust = 1) :
            isContinuation(false),
            g(_g),
            numThreads(nthr),
            process_rank(rank),
            size_of_cluster(nclust) {}

    virtual ~RootJohnTask() {}

    virtual tbb::task *execute() override;

protected:
    bool isContinuation;
    int numThreads;
    Graph *g;

    int process_rank;
    int size_of_cluster;
};

tbb::task *RootJohnTask::execute() {

    if (!isContinuation) {
        set_ref_count(1);
        parallel_for(size_t(0), size_t(g->getVertexNo()),  [&](size_t i) {
            if ((g->numNeighbors(i) != 0) && (g->numInEdges(i) != 0)) {
                parallel_for(size_t(g->offsArray[i]), size_t(g->offsArray[i + 1]), [&](size_t ind) {
                    int w = g->edgeArray[ind].vertex;
                    auto &tset = g->edgeArray[ind].tstamps;

                   parallel_for(size_t(0), size_t(tset.size()), [&](size_t j) {
                       if((ind+j) % size_of_cluster == process_rank) {
                           BlockedSet *blocked = new BlockedSet(g->getVertexNo());
                           BlockedList *BMap = new BlockedList;
                           Cycle *cycle = new Cycle();
                           cycle->push_back(i);

                           ThreadDataGuard *thrData = new ThreadDataGuard(cycle, blocked, BMap);

                           JohnsonsTask *a = new(allocate_child()) JohnsonsTask(g, w, tset[j], cycle, blocked, BMap,
                                                                                thrData);

                           increment_ref_count();
                           tbb::task::spawn(*a);
                       }
                    });
                });
            }
        });

        recycle_as_safe_continuation();
        isContinuation = true;
    }
    return NULL;
}
}

/// ************ Johnson algorithm with time window - top level ************

void allCyclesJohnsonCoarseGrainedTW(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    tbb::combinable<CycleHist> pt_chist_jh;

    //int start = process_rank;
    parallel_for(size_t(0), size_t(g->getVertexNo()), [&](size_t i) {
        if ((g->numNeighbors(i) != 0) && (g->numInEdges(i) != 0)) {
            parallel_for(size_t(g->offsArray[i]), size_t(g->offsArray[i + 1]), [&](size_t ind) {

                int w = g->edgeArray[ind].vertex;
                auto &tset = g->edgeArray[ind].tstamps;

                parallel_for(size_t(0), size_t(tset.size()), [&](size_t j) {

                    if((ind+j) % size_of_cluster == process_rank) {
                        auto &my_hist = pt_chist_jh.local();

                        BlockedSet blocked(g->getVertexNo());
                        BlockedList BMap;
                        Cycle *cycle = new Cycle();
                        cycle->push_back(i);

                        cyclesJohnsonTW(g, w, tset[j], cycle, blocked, BMap, my_hist);

                        delete cycle;
                    }
                });

            });
        }
    });

    combineCycleHistogram(pt_chist_jh, result);
}

void allCyclesJohnsonFineGrainedTW(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    tbb::task *a = new(RootJohnTask::allocate_root()) RootJohnTask(g, numThreads, process_rank, size_of_cluster);
    tbb::task::spawn_root_and_wait(*a);

    combineCycleHistogram(pt_chist_jh, result);
}