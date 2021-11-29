#ifndef _CYCLE_ENUMERATION_
#define _CYCLE_ENUMERATION_

#include <map>
#include <list>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "Graph.h"
#include "DataStructs.h"
#include "Macros.h"

#include <tbb/task.h>
#include <tbb/atomic.h>
#include <tbb/combinable.h>
#include <tbb/mutex.h>
#include <tbb/tick_count.h>
#include <tbb/spin_mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

using namespace std;

struct TstampInt {
    int first_ts;
    TimestampSet::iterator it_begin, it_end; // Indices of starting and ending timestamp
    TstampInt(int f, TimestampSet::iterator  b, TimestampSet::iterator  e)
        : first_ts(f), it_begin(b), it_end(e) {}
};

typedef ConcurrentList<int> Cycle;
typedef ConcurrentList<TstampInt> TimestampGroups;
typedef HashSet BlockedSet;
typedef unordered_map<int, HashSet> BlockedList;
typedef map<int, uint64_t> CycleHist;
typedef HashSet StrongComponent;
typedef list<Cycle*> CycleList;

/// Recording cycles
void recordCycle(Cycle* current, CycleHist& result);
void recordBundledCycle(Cycle* current, TimestampGroups* tg, CycleHist& result);
void combineCycleHistogram(tbb::combinable<CycleHist>& pt_chist, CycleHist& result);

struct TempEdge {
    int vertex;
    TimestampSet::iterator ts_begin, ts_end; // Iterators of starting and ending timestamp
    TempEdge() : vertex(-1) {}
    TempEdge(int v, TimestampSet::iterator t0, TimestampSet::iterator t1)
        : vertex(v), ts_begin(t0), ts_end(t1) {}
};

bool edgeInTimeInterval(int tstart, int timeWindow, int vstart,  int vert, TimestampSet& tset);

typedef VectorPath<TempEdge> TempPathBundle;
typedef VectorPath<int> Path;
typedef HashMap BlockedMap;

/// Enumerating simple cycles within time window
void allCyclesJohnsonCoarseGrainedTW(Graph *g, CycleHist& result, int numThreads);
void allCyclesJohnsonFineGrainedTW(Graph *g, CycleHist& result, int numThreads);

void allCyclesReadTarjanFineGrainedTW(Graph *g, CycleHist& result, int numThreads);
void allCyclesReadTarjanCoarseGrainedTW(Graph *g, CycleHist& result, int numThreads);

/// Enumerating temporal cycles
void allCyclesReadTarjanTemp(Graph *g, CycleHist& result);
void allCyclesTempReadTarjanCoarseGrained(Graph *g, CycleHist& result, int numThreads);
void allCyclesTempReadTarjanFineGrained(Graph *g, CycleHist& result, int numThreads);

void allCyclesTempJohnsonFineGrained(Graph *g, CycleHist& result, int numThreads);
void allCyclesTempJohnsonCoarseGrained(Graph *g, CycleHist& result, int numThreads);

/// Preprocessing functions
int findCycleUnions(Graph *g, EdgeData startEdge, int startVert, int timeWindow, StrongComponent *&cunion);
void cycleUnionExecTime(Graph *g, int numThreads);

/// dfs subroutine
bool findPathTemp(Graph *g, TempEdge e, EdgeData start, HashMapStack& blocked, TempPathBundle*& path, StrongComponent *cunion = NULL, bool nopath = false);
bool dfsPruneTemp(Graph *g, TempEdge e, EdgeData start, HashMapStack& blocked, TempPathBundle*& path, StrongComponent *cunion = NULL, bool nopath = false);
bool findPath(Graph *g, int u, int start, HashSetStack& blocked, Path*& path, int tstart = -1);
bool dfsPrune(Graph *g, int u, int start, HashSetStack& blocked, Path*& path, int tstart = -1);

/// 2SCENT
struct TmpEdge {
    int from, to, ts;
    TmpEdge(int _f = 0, int _t = 0, int _ts = 0) : from(_f), to(_t), ts(_ts) {}
};

typedef unordered_map<int, int> UnblockList;
typedef unordered_map<int, UnblockList> UnblockLists;
typedef unordered_map<int, unordered_set<int>> SummarySet;
typedef unordered_map<int, SummarySet> SummarySets;
typedef vector<TmpEdge> EdgeList;

class ClosingTimes {
private:
    unordered_map<int, int> elems;
public:
    ClosingTimes(int size = 0) {}
    void insert(int node, int ts) {
        if(ts == std::numeric_limits<int>::max()) elems.erase(node);
        else elems[node] = ts;
    }
    int at(int node) {
        auto it = elems.find(node);
        if(it == elems.end()) return std::numeric_limits<int>::max();
        return it->second;
    }
    void clear() { elems.clear();}
    int size() { return elems.size(); }
    unordered_map<int, int>::iterator begin() { return elems.begin(); }
    unordered_map<int, int>::iterator end() { return elems.end(); }
};

struct Seed {
    int root, ts, te, tn;
    unordered_set<int> cands;
    Seed(int _r = -1, int _ts = -1, int _te = -1) : root(_r), ts(_ts), te(_te), tn(-1) {}
};

typedef unordered_map<int, vector<Seed>> SeedSet;

void allCycles2scent(Graph *g, CycleHist& result);
void allCycles2scentSeed(Graph *g, SeedSet& seeds, CycleHist& result);
double createSeeds2scent(Graph *g, SeedSet*& result);
void seedTest2scent(Graph *g);


#endif