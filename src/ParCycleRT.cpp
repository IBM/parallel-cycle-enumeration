#include "CycleEnumeration.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <vector>
#include <chrono>
#include <tbb/task.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/mutex.h>

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;
using namespace tbb;

tbb::combinable<CycleHist> pt_rtchist;

extern int timeWindow;

//********************* Main task ****************************

namespace {

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
            blocked = new HashSetStack(true /*concurrent*/);
            #endif
            current = new Cycle(true /*concurrent*/);
        }

        ThreadDataGuard(ThreadDataGuard* guard, int lvl, int pathSize) : refCount(1) {
            #ifdef BLK_FORWARD
            blocked = guard->blocked->clone(lvl);
            #endif
            current = guard->current->clone(pathSize);
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
                delete current; current = NULL;
                delete this;
            }
        }

        #ifdef BLK_FORWARD
        HashSetStack* blocked = NULL;
        #endif
        Cycle* current = NULL;

        int refCount = 1;

        typedef tbb::spin_mutex RefCountMutexType;
        RefCountMutexType RefCountMutex;
    };

    class RTCycleTask : public tbb::task {
    public:
        RTCycleTask(Graph *_g, EdgeData _e, Cycle *_current, HashSetStack *_blocked,
                    ThreadDataGuard *tdg, Path *main, Path *alt, bool root = false, int lvl = 0, int tst = -1) :
                g(_g), e(_e), current(_current), isRoot(root), level(lvl), thrData(tdg), tstart(tst),
                blocked(_blocked),  pathSize(_current->size()), ind(g->offsArray[e.vertex]),
                main_path(main), alt_path(alt) {}

        virtual ~RTCycleTask() {}

        virtual tbb::task *execute() override;

//protected:

        bool isContinuation = false;

        tbb::task *Continuation();

        // sequential execution
        tbb::task *CycleReadTarjan();

        //Parameters
        Graph *g = NULL; EdgeData e; Cycle *current = NULL;
        int pathSize, tstart;
        HashSetStack *blocked = NULL;
        Path *main_path, *alt_path;
        bool isRoot = false;
        ThreadDataGuard *thrData;

        /// Task control
        int level = 0;
        bool rootTask = false;
        int ownerThread = -1;

        /// Continuation stealing
        bool newTask = true;
        int ind = 0;

#ifdef NEIGHBOR_ORDER
        vector<Path*> allPaths;
#endif
        int it = -1;
    };

    class RTFollowPathTask : public tbb::task {
    public:
        RTFollowPathTask(Graph *_g, EdgeData _e, Cycle *_current, HashSetStack *_blocked,
                         ThreadDataGuard *tdg, Path *main, bool root = false, int lvl = 0, int tst = -1) :
                g(_g), e(_e), current(_current), isRoot(root), thrData(tdg),
                blocked(_blocked), level(lvl), tstart(tst),
                current_path(main) {}

        virtual ~RTFollowPathTask() {}

        virtual tbb::task *execute() override;

//protected:

        bool isContinuation = false;

        tbb::task *Continuation();

        // sequential execution
        tbb::task *FollowPath();

        ThreadDataGuard *thrData;

        //Parameters
        Graph *g = NULL; EdgeData e; int tstart;
        Cycle *current = NULL;
        HashSetStack *blocked = NULL;
        Path *current_path = NULL;

        bool isRoot = false;
        int level = 0;
    };

/************* execute *******************/
    tbb::task *RTCycleTask::execute() {
        if(newTask) { set_ref_count(1); newTask = false;}

        if (!isContinuation) {
            int thisThreadId = getTheadId();
            if(ownerThread != thisThreadId) {
                ownerThread = thisThreadId;

                // Copy on steal
                ThreadDataGuard* newb = new ThreadDataGuard(thrData, level, pathSize);

                // Decrement the ref. count of the previous blocked map
                thrData->decrementRefCount();

                thrData = newb;
                // Update the pointers
                #ifdef BLK_FORWARD
                blocked = thrData->blocked;
                #endif
                current = thrData->current;
            }
            else
            {
                #ifdef BLK_FORWARD
                blocked->setLevel(level);
                #endif
                current->pop_back_until(pathSize);
            }

            tbb::task *nextTask = CycleReadTarjan();

            return nextTask;
        } else return Continuation();
    }

    tbb::task *RTFollowPathTask::execute() {
        if (!isContinuation) {
            set_ref_count(1);

            tbb::task *nextTask = FollowPath();

            recycle_as_safe_continuation();
            isContinuation = true;
            return nextTask;
        } else
            return Continuation();
    }

    tbb::task *RTCycleTask::CycleReadTarjan() {

        tbb::task* retTask = NULL;


        #if defined(NEIGHBOR_ORDER)
        int ind = g->offsArray[e.vertex];
        if(it == -1)
        #endif
        while(ind < g->offsArray[e.vertex+1]) {
            int w = g->edgeArray[ind].vertex;
            auto &tset = g->edgeArray[ind].tstamps;
            ind++;

            if (tstart != -1) {
                if (!edgeInTimeInterval(tstart, timeWindow, current->front(), e.vertex, tset)) continue;
            }

            if ((tstart == -1) && (w < current->front())) continue;
            if (blocked->exists(w)) continue;

            Path *current_path = NULL;
            if (main_path && w == main_path->back()) {
                current_path = main_path;
                main_path = NULL;
            } else if (alt_path && w == alt_path->back()) {
                current_path = alt_path;
                alt_path = NULL;
            } else {
                bool found = false;

                found = findPath(g, w, current->front(), *blocked, current_path, tstart);

                if (!found) continue;
            }


    #ifdef NEIGHBOR_ORDER
            allPaths.push_back(current_path);
    }

    if(it == -1) it = 0;
    while(it < allPaths.size()) {
        Path* current_path = allPaths[it];
        it++;
    #endif

            // Forwarding the blocked set
            #ifdef BLK_FORWARD
            HashSetStack* new_blocked = blocked;
            new_blocked->incrementLevel();
            #else
            HashSetStack *new_blocked = new HashSetStack(g->getVertexNo());
            for (auto c : *current) if (c != current->front()) new_blocked->insert(c);
            #endif

            RTFollowPathTask *a = new(allocate_child()) RTFollowPathTask(g, EdgeData(e.vertex, -1), current, new_blocked,
                                                                  thrData, current_path, false, level, tstart);

            increment_ref_count();
            #if !defined(NEIGHBOR_ORDER)
            bool exitCond = (ind == g->offsArray[e.vertex + 1]);
            #else
            bool exitCond = (it == allPaths.size());
            #endif
            if (!exitCond) {
                recycle_to_reexecute();
                thrData->incrementRefCount();
                return a;
            } else retTask = a;

        }


    if(!retTask) {
        thrData->decrementRefCount();
    }
    #ifndef BLK_FORWARD
    delete blocked; blocked = NULL;
    #endif

        recycle_as_safe_continuation();
        isContinuation = true;

        return retTask;
    }

    tbb::task *RTFollowPathTask::FollowPath() {
        tbb::task *a = NULL;

        Path *another_path = NULL;
        bool branching = false;

        int prev_vertex = -1;
        while (current_path->back() != current->front() && !branching) {
            prev_vertex = current_path->back();
            current_path->pop_back();
            current->push_back(prev_vertex);
            blocked->insert(prev_vertex);

            for (int ind = g->offsArray[prev_vertex]; ind < g->offsArray[prev_vertex + 1]; ind++) {
                int u = g->edgeArray[ind].vertex;
                auto& tset = g->edgeArray[ind].tstamps;

                if(tstart != -1) {
                    if(!edgeInTimeInterval(tstart, timeWindow, current->front(), prev_vertex, tset)) continue;
                }

                if (u != current_path->back() && ((tstart != -1) || (tstart == -1) && (u > current->front())) && !blocked->exists(u)
                    && !branching) {

                    branching = dfsPrune(g, u, current->front(), *blocked, another_path, tstart);

                    if (branching) {
                        #ifndef BLK_FORWARD
                        delete another_path; another_path = NULL;
                        #endif
                        break;
                    }
                }
            }
        }

        if (branching) {
            #ifndef BLK_FORWARD
            if (current_path) delete current_path; current_path = NULL;
            #endif
            increment_ref_count();
            a = new(allocate_child()) RTCycleTask(g, EdgeData(prev_vertex, -1), current, blocked, thrData, current_path,
                                                  another_path, false, level+1, tstart);
        } else {

            delete current_path; current_path = NULL;
            auto &my_hist = pt_rtchist.local();
            recordCycle(current, my_hist);

            #if !defined(BLK_FORWARD)
            delete blocked; blocked = NULL;
            #endif

            thrData->decrementRefCount();
        }
        return a;
    }

    tbb::task *RTCycleTask::Continuation() {
        return NULL;
    }

    tbb::task *RTFollowPathTask::Continuation() {
        return NULL;
    }
}
//********************* Root task ****************************

class RootRTCycleTask : public tbb::task {
public:
    RootRTCycleTask(Graph *g, int nthr, int rank = 0, int csize = 1):
            isContinuation(false), process_rank(rank), size_of_cluster(csize), graph(g), numThread(nthr){}

    void rootAllSimpleCyclesTW();

    virtual tbb::task* execute() override;

protected:
    bool isContinuation;
    Graph *graph;
    int numThread;
    int process_rank = 0;
    int size_of_cluster = 1;
};

tbb::task* RootRTCycleTask::execute() {
    if(!isContinuation) {
        set_ref_count(1);
        rootAllSimpleCyclesTW();

        recycle_as_safe_continuation();
        isContinuation = true;
    }

    return NULL;
}

void RootRTCycleTask::rootAllSimpleCyclesTW() {
    parallel_for(size_t(0), size_t(graph->getVertexNo()), [&](size_t i) {
        if ((graph->numNeighbors(i) != 0) && (graph->numInEdges(i) != 0)) {

            parallel_for(size_t(graph->offsArray[i]), size_t(graph->offsArray[i + 1]), [&](size_t ind) {
            int w = graph->edgeArray[ind].vertex;
            auto &tset = graph->edgeArray[ind].tstamps;

            parallel_for(size_t(0), size_t(tset.size()), [&](size_t j) {

                if((ind+j) % size_of_cluster == process_rank) {
                    ThreadDataGuard *thrData = nullptr;
                    thrData = new ThreadDataGuard();

                    #if defined(BLK_FORWARD)
                    HashSetStack *blocked = thrData->blocked;
                    #else
                    HashSetStack *blocked = new HashSetStack(graph->getVertexNo());
                    #endif

                    Cycle *cycle = thrData->current; cycle->push_back(i);

                    Path *current_path = NULL;
                    bool found = findPath(graph, w, i, *blocked, current_path, tset[j]);

                    if (found) {
                        #if  defined(BLK_FORWARD)
                        blocked->incrementLevel();
                        #endif

                        increment_ref_count();
                        spawn(*new(allocate_child()) RTFollowPathTask(graph, EdgeData(i, -1), cycle, blocked,
                                thrData, current_path, false, 0, tset[j]));
                    }
                }
                });
            });
        }
    });
}

/// ************ fine-grained Read-Tarjan algorithm with time window - top level ************

void allCyclesReadTarjanFineGrainedTW(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    tbb::task *rt = new(tbb::task::allocate_root()) RootRTCycleTask(g, numThreads, process_rank, size_of_cluster);
    tbb::task::spawn_root_and_wait(*rt);

    combineCycleHistogram(pt_rtchist, result);
}
