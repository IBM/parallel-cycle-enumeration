#include "CycleEnumeration.h"
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <utility>
#include <ctime>
#include <chrono>
#include <queue>
#include <iostream>
#include <fstream>

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;
using namespace tbb;

extern int timeWindow;
extern bool useCUnion;
extern ConcurrentCounter vertexVisits;

///************************** 2SCENT implementation ************************************

void unblock2scent(Graph* g, int node, int lastts, ClosingTimes& ctime, UnblockLists& Ulists) {
    if(lastts > ctime.at(node)) {
        ctime.insert(node, lastts);

        vector<pair<int, int>> tmpArray;
        tmpArray.reserve(Ulists[node].size());

        auto it = Ulists[node].begin();
        while(it != Ulists[node].end()) {
            if(it->second < lastts) {
                int w = it->first;
                it = Ulists[node].erase(it);

                int new_blocked_ts = -1;
                int ts_max = -1;

                bool node_found = false;
                for(int ind = g->offsArray[w]; ind < g->offsArray[w+1]; ind++) {
                    int u = g->edgeArray[ind].vertex;
                    auto& tset = g->edgeArray[ind].tstamps;

                    for(auto ts : tset) {
                        if (u == node) {
                            node_found = true;
                            if (ts < lastts) ts_max = max(ts, ts_max);
                            if (ts >= lastts) {
                                if (new_blocked_ts == -1) new_blocked_ts = ts;
                                else new_blocked_ts = min(new_blocked_ts, ts);
                            }
                        }
                        if (node_found && u != node) break;
                    }
                }

                if(new_blocked_ts != -1) tmpArray.push_back(std::make_pair(w, new_blocked_ts));
                unblock2scent(g, w, ts_max, ctime, Ulists);
            }
            else it++;
        }

        for(auto edge : tmpArray) Ulists[node].insert(edge);
    }
}

void extend2scent(Graph* g, int node, int ts, UnblockList& Ul) {
    auto it = Ul.find(node);
    if(it == Ul.end() || it != Ul.end() && it->second > ts) {
        Ul[node] = ts;
    }
}

int cycles2scentBundled(Graph* g, EdgeData e, EdgeData start, Cycle* current, TimestampGroups* tg, ClosingTimes& ctime,
                         UnblockLists&  Ulists, CycleHist& result, Seed *seed = NULL) {

    current->push_back(e.vertex);

    if(tg->size() != 0) {
        int tmin = tg->back().first_ts;
        ctime.insert(e.vertex, tmin);
    }

    vertexVisits.increment();
    int lastp = 0;

    /// Iterate through the neighbors of e.vertex
    for(int ind = g->offsArray[e.vertex]; ind < g->offsArray[e.vertex+1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;

        /// Determining the timestamps in the edge bundle
        int prev_min_ts = tg->back().first_ts;

        int end_tstamp = start.tstamp + timeWindow;
        if(seed != NULL && seed->te != -1) end_tstamp = seed->te;

        auto it_start = upper_bound(tset.begin(), tset.end(), prev_min_ts);

        if(it_start == tset.end() || *it_start > end_tstamp) continue;

        if((seed != NULL) && (seed->cands.find(w) == seed->cands.end())) continue;

        int ts = *(it_start);

        /// Progress recursivelly
        if (w == current->front()) {
            auto it_end = upper_bound(it_start, tset.end(), end_tstamp);
            tg->push_back(TstampInt(*(it_start), it_start, it_end));
            recordBundledCycle(current, tg, result);
            tg->pop_back();

            int tmax = *(it_end-1);
            lastp = ((tmax > lastp) ? tmax : lastp);
        } else {
            if(*(it_start) < ctime.at(w)) {
                auto it_end = upper_bound(it_start, tset.end(), end_tstamp);
                auto it_ctime = lower_bound(it_start, tset.end(), ctime.at(w));

                tg->push_back(TstampInt(*(it_start), it_start, it_ctime));
                int lastx = cycles2scentBundled(g, EdgeData(w, ts), start, current, tg, ctime, Ulists, result, seed);
                tg->pop_back();

                // Determining tmin
                int tmin = -1;
                int new_lastp = -1;

                auto it_lastx = lower_bound(it_start, it_end, lastx);
                if(it_lastx != it_start) new_lastp = *(it_lastx-1);
                if(it_lastx != it_end) tmin = *it_lastx;

                lastp = ((new_lastp > lastp) ? new_lastp : lastp);

                if (tmin != -1) {
                    extend2scent(g, e.vertex, tmin, Ulists[w]);
                }
            }
            else if(*(it_start) >= ctime.at(w)) {
                extend2scent(g, e.vertex, *(it_start), Ulists[w]);
            }
        }
    }

    if(lastp > 0) unblock2scent(g, e.vertex, lastp, ctime, Ulists);
    current->pop_back();

    return lastp;
}

/// ********************** Fine-grained temporal Johnson algorithm ***********************

namespace {
tbb::combinable<CycleHist> pt_cycleHist;

/// Thread Control block
tbb::atomic<int> globalId = 0;
tbb::enumerable_thread_specific<int> threadIds;

int getThreadId() {
    bool exists = true;
    auto& myid = threadIds.local(exists);
    if(!exists) myid = globalId++;
    return myid;
}

typedef ConcurrentList<int> OldLockStack;

typedef tbb::spin_mutex ThrDataMutexType;

struct ThreadDataGuard {
public:
    ThreadDataGuard(Graph *g, Cycle* c, TimestampGroups* t, ClosingTimes* ct, UnblockLists* ul) :
            graph(g), cycle(c), tg(t), oldCtimes(new OldLockStack),
            ctime(ct), ulists(ul), refCount(1) {}

    ThreadDataGuard(ThreadDataGuard* guard, int pathSize) : graph(guard->graph), refCount(1) {
        // Copy the data from another thread
       {
            ThrDataMutexType::scoped_lock lock(guard->ThrDataMutex);
            ctime = new ClosingTimes(*(guard->ctime));
            ulists = new UnblockLists(*(guard->ulists));
            cycle = new Cycle(*(guard->cycle));
            tg = guard->tg->clone(pathSize-1);
            oldCtimes = new OldLockStack(*(guard->oldCtimes));
        }
        // Remove the invalid vertices
        while(cycle->size() > pathSize) {
            int lastVertex = cycle->back(); cycle->pop_back();
            int lastCtime = oldCtimes->back(); oldCtimes->pop_back();
            unblock2scent(graph, lastVertex, lastCtime, *ctime, *ulists);
        }

        for(int i = 1; i < pathSize; i++) {
            int vert = cycle->at(i);
            int tmin = tg->at(i-1).first_ts;
            ctime->insert(vert, tmin);
        }
    }

    void incrementRefCount() {
        ThrDataMutexType::scoped_lock lock(RefCountMutex);
        refCount++;
    }

    void decrementRefCount() {
        ThrDataMutexType::scoped_lock lock(RefCountMutex);
        refCount--;
        if(refCount <= 0) {
            delete cycle; delete tg; delete ctime; delete ulists; delete oldCtimes;
            cycle = NULL; tg = NULL; ctime = NULL; ulists = NULL; oldCtimes = NULL;
            delete this;
        }
    }

    // Global graph data
    Graph *graph;

    // Guarded data
    Cycle* cycle = NULL;
    TimestampGroups* tg = NULL;
    ClosingTimes* ctime = NULL;
    UnblockLists* ulists = NULL;

    // Old ctime values
    OldLockStack* oldCtimes = NULL;

    ThrDataMutexType ThrDataMutex;
private:
    // Reference counter
    int refCount = 1;

    ThrDataMutexType RefCountMutex;
};

class TempJohnsonsTask : public tbb::task {
public:
    TempJohnsonsTask(Graph* g, int v, TstampInt e, EdgeData s,
            Cycle* cyc, TimestampGroups* t, ClosingTimes* ct, UnblockLists*  ul,
            ThreadDataGuard *tdg, StrongComponent* cu = NULL, TempJohnsonsTask* par = NULL):
            graph(g), vert(v), ebundle(e), start(s), parent(par),
            cycle(cyc), tg(t), ctime(ct), pathSize(cyc->size()),
            cunion(cu), ulists(ul), ownerThread(getThreadId()),
            thrData(tdg), ind(g->offsArray[vert])
    {}

    virtual ~TempJohnsonsTask(){}
    virtual tbb::task* execute() override;

    void copyOnSteal();

    /// Return lastx value from the recursive call on edge e
    void returnValue(int lastx);

    TimestampSet::iterator my_it_end;

    tbb::task* SpawnTask();
    tbb::task* Continuation();

    //Parameters
    int vert = -1;
    TstampInt ebundle;
    EdgeData start;
    Cycle* cycle = NULL;
    TimestampGroups* tg = NULL;
    ClosingTimes* ctime = NULL;
    UnblockLists* ulists = NULL;

    Graph *graph = NULL;
    StrongComponent *cunion = NULL;

    //Return
    ThrDataMutexType ReturnMutex;
    TempJohnsonsTask* parent = NULL;
    tbb::atomic<int> retLastp = false;

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

tbb::task* TempJohnsonsTask::execute() {
    if (!isContinuation) return SpawnTask();
    else return Continuation();
}

void TempJohnsonsTask::copyOnSteal() {
    /// Copy-on-steal
    int thisThreadId = getThreadId();

    //if (is_stolen_task()) {
    if(ownerThread != thisThreadId) {
        ownerThread = thisThreadId;

        stolenTask = true;

        // Copy on steal
        ThreadDataGuard* newThrData = new ThreadDataGuard(thrData, pathSize);

        // Decrement the ref. count of the previous blocked map
        thrData->decrementRefCount();

        thrData = newThrData;
        // Update the pointers
        cycle = thrData->cycle;
        tg = thrData->tg;
        ctime = thrData->ctime;
        ulists = thrData->ulists;
    }
}

tbb::task* TempJohnsonsTask::SpawnTask() {
    copyOnSteal();

    /// If we are executing the task for the first time
    if(newTask) {
        set_ref_count(1);
        newTask = false;
        vertexVisits.increment();
        {
            ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
            cycle->push_back(vert);
            tg->push_back(ebundle);
            pathSize++;

            if(tg->size() != 0) {
                int tmin = ebundle.first_ts;
                int prevCtime = ctime->at(vert);
                ctime->insert(vert, tmin);
                thrData->oldCtimes->push_back(prevCtime);
            }
        }
    }

    tbb::task *retTask = NULL;

    /// Iterate through the neighbors of vert
    while (ind < graph->offsArray[vert + 1]) {
        int w = graph->edgeArray[ind].vertex;
        auto& tset = graph->edgeArray[ind].tstamps;
        ind++;

        int prev_min_ts = ebundle.first_ts;
        int end_tstamp = start.tstamp + timeWindow;

        /// Determining the starting timestamp
        auto it_start = upper_bound(tset.begin(), tset.end(), prev_min_ts);

        /// No timestamp fulfills the condition
        if(it_start == tset.end() || *it_start > end_tstamp) continue;

        /// Skip if a vertex is not in the cycle union
        if((cunion != NULL) && (!cunion->exists(w))) continue;

        int ts = *(it_start);

        /// Progress recursivelly
        if (w == start.vertex) {
            auto it_end = upper_bound(it_start, tset.end(), end_tstamp);
            {
                ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
                tg->push_back(TstampInt(*(it_start), it_start, it_end));
            }
            auto& local_hist = pt_cycleHist.local();
            recordBundledCycle(cycle, tg, local_hist);
            {
                ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
                tg->pop_back();
            }

            int tmax = *(it_end-1);
            {
                ThrDataMutexType::scoped_lock lock(ReturnMutex);
                if(tmax > retLastp) retLastp = tmax;
            }
        } else {

            if(*(it_start) < ctime->at(w)) {
                /// Spawning child task
                increment_ref_count();

                auto it_end = upper_bound(it_start, tset.end(), end_tstamp);
                auto it_ctime = lower_bound(it_start, it_end, ctime->at(w));

                auto ebundle = TstampInt(*(it_start), it_start, it_ctime);
                TempJohnsonsTask *a = new(allocate_child())  TempJohnsonsTask(graph, w, ebundle, start,
                        cycle, tg, ctime, ulists, thrData, cunion, this);

                a->my_it_end = it_end;

                /// Continuation stealing
                thrData->incrementRefCount();
                if (ind < graph->offsArray[vert + 1]) {
                    recycle_to_reexecute();
                    return a;
                } else retTask = a;
            }
            else if(*(it_start) >= ctime->at(w)) {
                ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);
                extend2scent(graph, vert, *(it_start), (*ulists)[w]);
            }
        }
    }

    // No continuation happens if the task is stolen
    if(stolenTask) thrData->decrementRefCount();

    recycle_as_safe_continuation();
    isContinuation = true;

    return retTask;
}

void TempJohnsonsTask::returnValue(int lastx) {
    ThrDataMutexType::scoped_lock lock(ReturnMutex);
    if (lastx > retLastp) retLastp = lastx;
}

tbb::task* TempJohnsonsTask::Continuation() {
    if(!stolenTask) {
        if(getThreadId() == ownerThread && pathSize == cycle->size()) {
            /// Compute tmin and new_lastp
            int prevVert = cycle->at(cycle->size()-2);

            int tmin = -1;
            int new_lastp = -1;
            auto it_start = ebundle.it_begin;

            auto it_end = my_it_end;
            auto it_lastx = lower_bound(it_start, it_end, retLastp);

            if(it_lastx != it_start) new_lastp = *(it_lastx-1);
            if(it_lastx != it_end) tmin = *it_lastx;

            /// Atomic update
            {
                ThrDataMutexType::scoped_lock lock(thrData->ThrDataMutex);

                cycle->pop_back();
                pathSize--;
                thrData->oldCtimes->pop_back();
                tg->pop_back();

                if (retLastp > 0) {
                    unblock2scent(graph, vert, retLastp, *ctime, *ulists);
                }

                if (tmin != -1) extend2scent(graph, prevVert, tmin, (*ulists)[vert]);
            }

            /// Return
            if (parent) parent->returnValue(new_lastp);

        }

        thrData->decrementRefCount();
    }

    if(!parent && cunion) { delete cunion; cunion = NULL; }
    return NULL;
}
}

/// ********************* creating seeds - 2SCENT *********************

bool sortFun (TmpEdge a, TmpEdge b) { return (a.ts < b.ts); }
bool seedSortFun (Seed a, Seed b) { if(a.ts == b.ts) return (a.te < b.te); else return (a.ts < b.ts); }

EdgeList createEdgeList(Graph *g) {
    EdgeList elist;
    elist.reserve(g->tstampNo);
    for (int from = 0; from < g->getVertexNo(); from++) {

        for (int ind = g->offsArray[from]; ind < g->offsArray[from + 1]; ind++) {
            int to = g->edgeArray[ind].vertex;
            auto& tset = g->edgeArray[ind].tstamps;
            for(auto ts : tset) {
                elist.push_back(TmpEdge(from, to, ts));
            }
        }
    }

    std::sort(elist.begin(), elist.end(), sortFun);

    return elist;
}

SeedSet generateSeeds(EdgeList& elist, int tw = timeWindow) {
    SummarySets Ssets;
    SeedSet result;

    for(auto edge : elist) {
        int a = edge.from, b = edge.to;
        SummarySet& SetA = Ssets[a];
        SummarySet& SetB = Ssets[b];

        SetB[a].insert(edge.ts);

        if(SetA.size() != 0) {
            auto it = SetA.begin();
            while(it != SetA.end()) {

                int x = it->first;
                auto& tstampSet = it->second;

                auto its = tstampSet.begin();
                while (its != tstampSet.end()) {
                    int tx = *its;
                    if (tx < edge.ts - tw) its = tstampSet.erase(its);
                    else its++;
                }

                if (tstampSet.size() == 0) it = SetA.erase(it);
                else it++;
            }

            for(auto sset : SetA) {
                int x = sset.first;
                for(auto tx : sset.second) SetB[x].insert(tx);
            }

            for(auto tb : SetB[b]) {
                Seed seed(b, tb, edge.ts);
                seed.tn = seed.ts + tw;
                seed.cands.insert(b);

                for(auto tedge : SetB) {
                    int c = tedge.first;
                    for(auto tc : tedge.second) {
                        if(tc > tb) {
                            seed.cands.insert(c);
                            break;
                        }
                    }
                }

                result[b].push_back(seed);
            }

            auto itb = SetB.find(b);
            SetB.erase(itb);
        }
    }

    return result;
}

SeedSet* combineSeeds(SeedSet& seeds, int tw = timeWindow) {
    SeedSet* compressed = new SeedSet;

    for(auto& sset : seeds) {
        vector<Seed>& seedVec = sset.second;

        Seed current = seedVec[0];
        int ts_cur = current.ts, te_cur = current.te;

        if(seedVec.size() == 1) { (*compressed)[sset.first].push_back(current); continue; }

        for(int i = 1; i < seedVec.size(); i++) {
            te_cur = seedVec[i].te;
            if(te_cur < ts_cur + tw) {
                if(current.te < te_cur) current.te = te_cur;
                current.cands.insert(seedVec[i].cands.begin(), seedVec[i].cands.end());
            } else {
                current.tn = seedVec[i].ts;
                (*compressed)[sset.first].push_back(current);

                current = seedVec[i];
                ts_cur = current.ts;
                te_cur = current.te;
            }

            if(i+1 == seedVec.size()) (*compressed)[sset.first].push_back(current);
        }
    }

    return compressed;
}

double createSeeds2scent(Graph *g, SeedSet*& result) {
    EdgeList elist = createEdgeList(g);

    auto seed_start = chrono::steady_clock::now();
    SeedSet seeds = generateSeeds(elist);
    for(auto& root : seeds) std::sort(root.second.begin(), root.second.end(), seedSortFun);
    result =  combineSeeds(seeds, timeWindow);
    auto seed_end = chrono::steady_clock::now();
    double seed_total = chrono::duration_cast<chrono::milliseconds>(seed_end - seed_start).count()/1000.0;

    return seed_total;
}

/// ************************** 2SCENT with and without preprocessing - top level **************************

void allCycles2scentSeed(Graph *g, SeedSet& seeds, CycleHist& result) {
    for(auto& seedVec : seeds) {

        int startVec = seedVec.first;

        for(auto& seed : seedVec.second) {

            int tstart = seed.ts, tend = seed.te, tn = seed.tn;
            unordered_set<int>& candVert = seed.cands;

            Cycle *current = new Cycle();
            current->push_back(startVec);

            ClosingTimes ctime(g->getVertexNo());
            UnblockLists Ulists;

            bool skipVertex = false;
            int prev_w = -1;

            for (int ind = g->offsArray[startVec]; ind < g->offsArray[startVec + 1]; ind++) {
                int w = g->edgeArray[ind].vertex;
                auto& tset = g->edgeArray[ind].tstamps;

                if (candVert.find(w) == candVert.end()) continue;

                auto it_start = lower_bound(tset.begin(), tset.end(), tstart);
                auto it_end = lower_bound(tset.begin(), tset.end(), tn);

                if(it_start >= it_end) continue;

                int ts = *it_start;

                TimestampGroups *tg = new TimestampGroups();
                tg->push_back(TstampInt(ts, it_start, it_end));
                cycles2scentBundled(g, EdgeData(w, ts), EdgeData(startVec, ts), current, tg, ctime, Ulists, result, &seed);
                delete tg;
            }

            delete current;
        }
    }
}

void allCycles2scent(Graph *g, CycleHist& result) {
    for (int i = 0; i < g->getVertexNo(); i++) {
        if(!(g->numNeighbors(i) == 0) && !(g->numInEdges(i) == 0)) {

            for(int ind = g->offsArray[i]; ind < g->offsArray[i+1]; ind++) {
                int w = g->edgeArray[ind].vertex;
                auto &tset = g->edgeArray[ind].tstamps;

                for (auto it = tset.begin(); it != tset.end(); it++) {
                    int ts = *it;

                    StrongComponent *cunion = NULL;
                    Seed *seed = NULL;

                    if (useCUnion) {
                        findCycleUnions(g, EdgeData(w, ts), i, timeWindow, cunion);
                        seed = new Seed;
                        for (auto el : *cunion) seed->cands.insert(el);
                    }

                    if (!cunion || cunion->size() > 1) {
                        ClosingTimes ctime(g->getVertexNo());
                        UnblockLists Ulists;
                        Cycle *current = new Cycle();
                        current->push_back(i);

                        TimestampGroups *tg = new TimestampGroups();
                        tg->push_back(TstampInt(ts, it, it+1));
                        int lastx = cycles2scentBundled(g, EdgeData(w, ts), EdgeData(i, ts), current, tg, ctime, Ulists, result, seed);
                        delete tg;

                        delete current;
                    }

                    if (cunion) {
                        delete cunion; delete seed;
                    }
                }
            }
        }
    }

}

/// ************************** Parallel Temporal Johnson - top level **************************

namespace {
class OuterLoopTempJohnTask : public tbb::task {
public:
    OuterLoopTempJohnTask(Graph *_g, tbb::combinable <CycleHist> &_ch, int _vert, int _w, int _ts,
                          TimestampSet::iterator it, bool tpar = false) :
            g(_g), pt_chist(_ch), vert(_vert), taskPar(tpar),
            w(_w), ts(_ts), iter(it) {}

    virtual ~OuterLoopTempJohnTask() {}

    tbb::task*  runCoarseGrained();
    tbb::task*  runFineGrained();

    virtual tbb::task *execute() override;

protected:
    int vert, w, ts;
    tbb::combinable <CycleHist> &pt_chist;
    TimestampSet::iterator iter;
    Graph *g;

    bool isContinuation = false;
    bool taskPar = false;
};

tbb::task* OuterLoopTempJohnTask::execute() {
    tbb::task* next = NULL;
    if(!isContinuation) {
        set_ref_count(1);

        if (!taskPar) next = runCoarseGrained();
        else next = runFineGrained();
    }

    return next;
}

tbb::task* OuterLoopTempJohnTask::runCoarseGrained() {
    auto &my_hist = pt_chist.local();

    StrongComponent *cunion = NULL;
    Seed *seed = NULL;

    //Filter-out vertices not belonging to the same cycle-union
    if (useCUnion) {
        findCycleUnions(g, EdgeData(w, ts), vert, timeWindow, cunion);
        seed = new Seed;
        for (auto el : *cunion) seed->cands.insert(el);
    }

    if (!cunion || cunion->size() > 1) {
        ClosingTimes ctime(g->getVertexNo());
        UnblockLists Ulists;
        Cycle *current = new Cycle();
        current->push_back(vert);

        TimestampGroups *tg = new TimestampGroups();
        tg->push_back(TstampInt(ts, iter, iter + 1));
        int lastx = cycles2scentBundled(g, EdgeData(w, ts), EdgeData(vert, ts), current, tg, ctime, Ulists, my_hist, seed);
        delete tg;

        delete current;
    }

    if (seed) {
        delete seed;
        delete cunion;
    }

    return NULL;
}

tbb::task* OuterLoopTempJohnTask::runFineGrained() {

    StrongComponent *cunion = NULL;

    //Filter-out vertices not belonging to the same cycle-union
    if (useCUnion) {
        findCycleUnions(g, EdgeData(w, ts), vert, timeWindow, cunion);
        if(cunion->size() <= 1) {
            delete cunion;
            return NULL;
        }
    }

    ClosingTimes *ctime = new ClosingTimes;
    UnblockLists *ulists = new UnblockLists;
    TimestampGroups *tg = new TimestampGroups();
    Cycle *cycle = new Cycle;
    cycle->push_back(vert);

    ThreadDataGuard *thrData = new ThreadDataGuard(g, cycle, tg, ctime, ulists);

    TstampInt ebundle = TstampInt(ts, iter, iter + 1);

    increment_ref_count();
    tbb::task * a = new(allocate_child())
        TempJohnsonsTask(g, w, ebundle, EdgeData(vert, ts), cycle, tg, ctime, ulists, thrData, cunion, NULL);

    recycle_as_safe_continuation();
    isContinuation = true;

    return a;
}

class RootTempJohnTask : public tbb::task {
public:
    RootTempJohnTask(Graph *_g, int nthr, tbb::combinable <CycleHist> &_ch, bool tpar = false, int rank = 0, int nclust = 1) :
            isContinuation(false), taskPar(tpar),
            g(_g), pt_chist(_ch),
            numThreads(nthr),
            process_rank(rank),
            size_of_cluster(nclust) {}

    virtual ~RootTempJohnTask() {}

    virtual tbb::task *execute() override;

protected:
    bool isContinuation;
    bool taskPar = false;
    int numThreads;
    Graph *g;
    tbb::combinable <CycleHist> &pt_chist;

    int process_rank;
    int size_of_cluster;
};

tbb::task *RootTempJohnTask::execute() {

    if (!isContinuation) {
        set_ref_count(1);

        parallel_for(size_t(0), size_t(g->getVertexNo()), [&](size_t i) {
            // Filter-out vertices having 0 in or out degree
            if (!(g->numNeighbors(i) == 0) && !(g->numInEdges(i) == 0)) {

                parallel_for(size_t(g->offsArray[i]), size_t(g->offsArray[i + 1]), [&](size_t ind) {

                    int w = g->edgeArray[ind].vertex;
                    auto &tset = g->edgeArray[ind].tstamps;

                    parallel_for(size_t(0), size_t(tset.size()), [&](size_t j) {

                        if((ind+j) % size_of_cluster == process_rank) {
                            int ts = tset[j];

                            increment_ref_count();
                            tbb::task *a = new(allocate_child()) OuterLoopTempJohnTask(g, pt_chist, i, w, ts, tset.begin() + j, taskPar);
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

void allCyclesTempJohnsonCoarseGrained(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    tbb::combinable<CycleHist> pt_cycleHist;

    task_scheduler_init init(numThreads);
    tbb::task *rt = new(RootTempJohnTask::allocate_root()) RootTempJohnTask(g, numThreads, pt_cycleHist, false, process_rank, size_of_cluster);
    tbb::task::spawn_root_and_wait(*rt);

    combineCycleHistogram(pt_cycleHist, result);
}

void allCyclesTempJohnsonFineGrained(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    tbb::task *rt = new(RootTempJohnTask::allocate_root()) RootTempJohnTask(g, numThreads, pt_cycleHist, true, process_rank, size_of_cluster);
    tbb::task::spawn_root_and_wait(*rt);

    combineCycleHistogram(pt_cycleHist, result);
}