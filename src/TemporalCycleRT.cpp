#include "CycleEnumeration.h"
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <ctime>
#include <chrono>
#include <queue>
#include <list>

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;

extern int timeWindow;
extern bool useCUnion;
extern ConcurrentCounter vertexVisits;

/// ********************** DFS functions for pruning the search space **********************

namespace {
// Prune the search space and find the path using DFS
    bool dfsUtilTemp(Graph *g, EdgeData u, EdgeData start, int depth, HashMapStack& blocked, BlockedMap& visited,
                     bool& blk_flag, TempPathBundle*& path, unsigned long& stat_cnt, StrongComponent *cunion = NULL, bool nopath = false) {

        blk_flag = false;
        bool this_blk_flag = true;

        if(u.vertex == start.vertex && depth == 0) {
            if(!nopath) path = new TempPathBundle(1);
            return true;
        }

        stat_cnt++;
#if defined(SUCCESSFUL_DFS_BLK)
        int prev_close = visited.at(u.vertex);
#endif
        if(u.vertex != start.vertex) visited.insert(u.vertex, -1);

        int maxts = -1;
        // We assume visited.at(u.vertex) cannot be -1 here
        int end_tstamp = blocked.at(u.vertex) == 0 ? start.tstamp + timeWindow : blocked.at(u.vertex);

        for(int ind = g->offsArray[u.vertex]; ind < g->offsArray[u.vertex+1]; ind++) {
            int w = g->edgeArray[ind].vertex;
            auto& tset = g->edgeArray[ind].tstamps;

            auto it_start = upper_bound(tset.begin(), tset.end(), u.tstamp);

            /// Update maxts
#if !defined(ONE_DAY_CYCLE) && !defined(NON_TEMPORAL)
            if(it_start != tset.begin()) {
                maxts = max(maxts, *(it_start-1));
            }
#else
            if (w < start.vertex) continue;
        /// MISSING
#endif

            /// If the timestamp for the given time window exists, or the vertex is not in the cycle-union
            if(it_start == tset.end() || *it_start > end_tstamp || (cunion && !cunion->exists(w))) continue;

            int ts = *it_start;

            if (w == start.vertex) {
                if (!nopath) {
                    path = new TempPathBundle(depth + 2);
                    auto it_end = upper_bound(it_start, tset.end(), end_tstamp);
                    path->push_back(TempEdge(w, it_start, it_end));
                }
                return true;
            }
#if defined(SUCCESSFUL_DFS_BLK)
                else if(visited.exists(w, ts)) {
            this_blk_flag = false;
        }
        else if (!blocked.exists(w, ts)) {
#else
            else if (!blocked.exists(w, ts) && !visited.exists(w, ts)) {
#endif
                bool rec_blk_flag = false;
                if (dfsUtilTemp(g, EdgeData(w, ts), start, depth + 1, blocked, visited, rec_blk_flag, path,
                                stat_cnt, cunion, nopath)) {
                    if (!nopath) {
                        auto it_end = upper_bound(it_start, tset.end(), end_tstamp);
                        path->push_back(TempEdge(w, it_start, it_end));
                    }
                    return true;
                }
                this_blk_flag &= rec_blk_flag;
            }
        }

        if (u.vertex != start.vertex) {
#if defined(SUCCESSFUL_DFS_BLK)
            if(this_blk_flag) {
            blocked.insert(u.vertex, maxts);
            visited.insert(u.vertex, prev_close);
        } else
#endif
            visited.insert(u.vertex, maxts);
        }
        blk_flag = this_blk_flag;
        return false;
    }

}

// Different wrappers for the purpose of collecting statistics and easier profiling
bool findPathTemp(Graph *g, TempEdge e, EdgeData start, HashMapStack& blocked,
                  TempPathBundle*& path, StrongComponent *cunion, bool nopath) {

    BlockedMap visited(g->getVertexNo());
    bool blck_flag = false;
    unsigned long vertex_visits = 0;

    EdgeData u(e.vertex, *(e.ts_begin));
    bool found = dfsUtilTemp(g, u, start, 0, blocked, visited, blck_flag, path, vertex_visits, cunion, nopath);
    if(found && !nopath) path->push_back(e);

    vertexVisits.increment(vertex_visits);

#ifdef BLK_FORWARD
    if (!found) blocked.include(visited);
#endif
    return found;
}

bool dfsPruneTemp(Graph *g, TempEdge e, EdgeData start, HashMapStack& blocked,
                  TempPathBundle*& path, StrongComponent *cunion, bool nopath)  {

    BlockedMap visited(g->getVertexNo());
    bool blck_flag = false;
    unsigned long vertex_visits = 0;

    EdgeData u(e.vertex, *(e.ts_begin));
    bool found = dfsUtilTemp(g, u, start, 0, blocked, visited, blck_flag, path, vertex_visits, cunion, nopath);
    if(found && !nopath) path->push_back(e);

    vertexVisits.increment(vertex_visits);

    if (!found) blocked.include(visited);
    return found;
}

/// ********************** Backtracking **********************
namespace {
void cyclesReadTarjanTemp(Graph* g, int vert, EdgeData start, Cycle* cycle, TimestampGroups* tg, int prevClose,
     HashMapStack& blocked, CycleHist& result, pair<TempPathBundle*, TempPathBundle*> paths, StrongComponent *cunion = NULL);

void followPathRT(Graph* g, int vert, EdgeData start, Cycle* cycle, TimestampGroups* tg, HashMapStack& blocked,
        CycleHist& result, TempPathBundle* current_path, StrongComponent *cunion = NULL) {
    TempPathBundle* another_path = NULL;
    bool branching = false;

    int prev_close = -1;
    int startVert = cycle->front();

    TempEdge prev;
    // Add vertices from the found path to the current path until a branching is found
    while(current_path->back().vertex != startVert && !branching) {
        prev = current_path->back();
        current_path->pop_back();

        /// Update the current cycle and path set
        cycle->push_back(prev.vertex);
        prev_close = blocked.at(prev.vertex);
        blocked.insert(prev.vertex, -1);

        int prev_min_ts = *(prev.ts_begin);
        int end_tstamp = prev_close == 0 ? start.tstamp + timeWindow : prev_close;

        /// Extract timestamps and add them to tg
        tg->push_back(TstampInt(prev_min_ts, prev.ts_begin, prev.ts_end));
        
        /// Detect branching and find an alternate path
        // For each neighboring vertex u
        for(int ind = g->offsArray[prev.vertex]; ind < g->offsArray[prev.vertex+1] && !branching; ind++) {
            int u = g->edgeArray[ind].vertex;
            auto& tset = g->edgeArray[ind].tstamps;

            #if defined(ONE_DAY_CYCLE) || defined(NON_TEMPORAL)
            if (u < startVert) continue;
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

            if ((!cunion || cunion->exists(u)) && !blocked.exists(u, tu) && !branching) {
                auto it_end = upper_bound(it_start, tset.end(), end_tstamp);

                #ifndef PATH_FORWARD
                branching = dfsPruneTemp(g, TempEdge(u, it_start, it_end), start, blocked, another_path, cunion, true);
                #else
                branching = dfsPruneTemp(g, TempEdge(u, it_start, it_end), start, blocked, another_path, cunion, false);
                #endif

                if (branching) break;
            }
        }
    }

    /// Record the cycles if needed
    if(!branching) {
        // Add the final timestamp group
        TempEdge finalEdge = current_path->back();
        tg->push_back(TstampInt(*(finalEdge.ts_begin), finalEdge.ts_begin, finalEdge.ts_end));

        delete current_path; current_path = NULL;
        recordBundledCycle(cycle, tg, result);

        tg->pop_back();
    }
        /// Recursive call
    else {
        #ifndef PATH_FORWARD
        if(current_path) delete current_path; current_path = NULL;
        #endif
        cyclesReadTarjanTemp(g, prev.vertex, start, cycle, tg, prev_close,
                             blocked, result, make_pair(current_path, another_path), cunion);

    }

    /// Remove all the vertices added after vert
    while (1) {
        int back = cycle->back();
        if(back == vert) break;
        cycle->pop_back();
        tg->pop_back();
    }

}

/// The root backtracking function
void cyclesReadTarjanTempRoot(Graph* g, int startVert, Cycle* cycle, TimestampGroups* tg, CycleHist& result) {
    /// Iterate through every neighboring vertex and find paths
    for(int ind = g->offsArray[startVert]; ind < g->offsArray[startVert+1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;

        for(auto it = tset.begin(); it != tset.end(); it++) {
            int tw = *it;

            #if defined(ONE_DAY_CYCLE) || defined(NON_TEMPORAL)
            if(w < startVert) continue;
            #endif

            EdgeData start(startVert, tw);

            StrongComponent *cunion = NULL;

            if (useCUnion) {
                findCycleUnions(g, EdgeData(w, tw), startVert, timeWindow, cunion);
            }

            if (!cunion || cunion->size() > 1) {
                /// Starting vertex and timestamp

                /// Find path that needs to be followed
                TempPathBundle *current_path = NULL;
                HashMapStack blocked;

                bool found = findPathTemp(g, TempEdge(w, it, it+1), start, blocked, current_path, cunion);
                if (!found) continue;

                TempPathBundle *another_path = NULL;
                bool branching = false;

                /// Traverse the path until we find a branching point
                followPathRT(g, startVert, start, cycle, tg, blocked, result, current_path, cunion);
            }
            delete cunion;
        }
    }
}

/// The main backtracking function
void cyclesReadTarjanTemp(Graph* g, int vert, EdgeData start, Cycle* cycle, TimestampGroups* tg, int prevClose,
     HashMapStack& blocked, CycleHist& result, pair<TempPathBundle*, TempPathBundle*> paths, StrongComponent *cunion) {
    int startVert = cycle->front();

#ifdef NEIGHBOR_ORDER
    vector<TempPathBundle*> allPaths;
#endif

    /// Iterate through every neighboring vertex and find paths
    for(int ind = g->offsArray[vert]; ind < g->offsArray[vert+1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;

        /// Find the next timestamp that will be used
        // Move index to the starting timestamp where tw >= tg->back().first_ts
        int prev_min_ts = tg->back().first_ts;
        int end_tstamp = prevClose == 0 ? start.tstamp + timeWindow : prevClose;

        auto it_start = upper_bound(tset.begin(), tset.end(), prev_min_ts);
        auto it_end = upper_bound(it_start, tset.end(), end_tstamp);

        if(it_start >= it_end) continue;

        int tw = *(it_start);

        /// Don't visit this vertex at all
        if ((cunion && !cunion->exists(w)) || blocked.exists(w, tw)) continue;

        /// Find path that needs to be followed
        TempPathBundle *current_path = NULL;
        if (paths.first && w == paths.first->back().vertex) {
            current_path = paths.first;
            paths.first = NULL;
        } else if (paths.second && w == paths.second->back().vertex) {
            current_path = paths.second;
            paths.second = NULL;
        } else {
            bool found = findPathTemp(g, TempEdge(w, it_start, it_end), start, blocked, current_path, cunion);
            if (!found) continue;
        }

#ifdef NEIGHBOR_ORDER
        allPaths.push_back(current_path);
    }

    for(auto current_path : allPaths) {
#endif
        /// Traverse the path until we find a branching point
#ifdef BLK_FORWARD
        blocked.incrementLevel();
        followPathRT(g, vert, start, cycle, tg, blocked, result, current_path, cunion);
        blocked.decrementLevel();
#else
        HashMapStack new_blocked;
        for (auto c : *cycle) if (c != startVert) new_blocked.insert(c, -1);
        followPathRT(g, vert, start, cycle, tg, new_blocked, result, current_path, cunion);
#endif
    }
    }
}

/// Top level function for root calls

void allCyclesReadTarjanTemp(Graph *g, CycleHist& result) {
    for(int i = 0; i < g->getVertexNo(); i++) {
        /// Filter-out vertices having 0 in or out degree
        if((g->numNeighbors(i) != 0) && (g->numInEdges(i) != 0)) {
            TimestampGroups *tsgroup = new TimestampGroups();
            Cycle *cycle = new Cycle();

            cycle->push_back(i);

            cyclesReadTarjanTempRoot(g, i, cycle, tsgroup, result);

            delete cycle;
            delete tsgroup;
        }
    }
}

void allCyclesTempReadTarjanCoarseGrained(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);

    tbb::combinable<CycleHist> pt_chist_rt;

    parallel_for(size_t(0), size_t(g->getVertexNo()), [&](size_t i) {
        /// Filter-out vertices having 0 in or out degree
        if(!(g->numNeighbors(i) == 0) && !(g->numInEdges(i) == 0)) {

            int startVert = i;
            /// Iterate through every neighboring vertex and find paths
            parallel_for(size_t(g->offsArray[startVert]), size_t(g->offsArray[startVert + 1]), [&](size_t ind) {
                int w = g->edgeArray[ind].vertex;
                auto& tset = g->edgeArray[ind].tstamps;

                parallel_for(size_t(0), size_t(tset.size()), [&](size_t it) {
                    int currRank = (ind + it) % size_of_cluster;
                    if (currRank == process_rank) {
                        int tw = tset[it];

                        StrongComponent *cunion = NULL;

                        if (useCUnion) findCycleUnions(g, EdgeData(w, tw), startVert, timeWindow, cunion);

                        if (!cunion || cunion->size() > 1) {
                            /// Starting vertex and timestamp
                            EdgeData start(startVert, tw);
                            TimestampGroups *tg = new TimestampGroups();
                            Cycle *cycle = new Cycle();
                            cycle->push_back(startVert);

                            /// Find path that needs to be followed
                            TempPathBundle *current_path = NULL;
                            HashMapStack blocked;

                            bool found = findPathTemp(g, TempEdge(w, tset.begin() + it, tset.begin() + it + 1), start, blocked, current_path, cunion);

                            if (found) {
                                /// Traverse the path until we find a branching point
                                auto &my_hist = pt_chist_rt.local();
                                followPathRT(g, startVert, start, cycle, tg, blocked, my_hist, current_path, cunion);
                            }

                            delete cycle; delete tg;
                        }
                        if (cunion) delete cunion;
                    }
                });
            });

        }
    });

    combineCycleHistogram(pt_chist_rt, result);
}
