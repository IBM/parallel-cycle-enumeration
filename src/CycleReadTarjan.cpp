#include "CycleEnumeration.h"
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <ctime>
#include <chrono>
#include <queue>
#include <list>

#ifdef MPI_IMPL
#include <mpi.h>
#endif

using namespace std;

extern int timeWindow;
extern ConcurrentCounter vertexVisits;

/// DFS/BFS functions for pruning the search space

// Prune the search space and find the path using DFS
bool dfsUtil(Graph *g, int u, int start, int depth, HashSetStack& blocked, BlockedSet& visited,
            bool& blk_flag, Path*& path, unsigned long& stat_cnt, int tstart = -1) {

    blk_flag = false;
    bool this_blk_flag = true;

    if(u == start && depth == 0) {
        path = new Path(1);
        path->push_back(u);
        return true;
    }

    stat_cnt++;
    visited.insert(u);

    for(int ind = g->offsArray[u]; ind < g->offsArray[u+1]; ind++) {
        int w = g->edgeArray[ind].vertex;
        auto& tset = g->edgeArray[ind].tstamps;

        if(tstart != -1) {
            if(!edgeInTimeInterval(tstart, timeWindow, start, u, tset)) continue;
        }

        if(w == start) {
            path = new Path(depth+2);
            path->push_back(w);
            path->push_back(u);
            return true;
        }
        else if(visited.exists(w)) {
            this_blk_flag = false;
        }
        else if(((tstart != -1) || (tstart == -1) && (w > start)) &&
            !blocked.exists(w) && !visited.exists(w)) {
            bool rec_blk_flag = false;
            if(dfsUtil(g, w, start, depth+1, blocked, visited, rec_blk_flag, path, stat_cnt, tstart)) {
                path->push_back(u);
                return true;
            }
            this_blk_flag &= rec_blk_flag;
        }
    }

#if defined(BLK_FORWARD) && defined(SUCCESSFUL_DFS_BLK)
    if (u != start) {
        if(this_blk_flag) {
            blocked.insert(u);
            visited.remove(u);
        }
    }
#endif
    blk_flag = this_blk_flag;
    return false;
}

// Different wrappers for the purpose of collecting statistics and easier profiling
bool findPath(Graph *g, int u, int start, HashSetStack& blocked, Path*& path, int tstart) {

    BlockedSet visited(g->getVertexNo());
    bool blck_flag = false;
    unsigned long vertex_visits = 0;
    bool found = dfsUtil(g, u, start, 0, blocked, visited, blck_flag, path, vertex_visits, tstart);

    vertexVisits.increment(vertex_visits);

#if defined(BLK_FORWARD) && defined(SUCCESSFUL_DFS_BLK)
    if (u != start && blck_flag) blocked.insert(u);
#endif

#ifdef BLK_FORWARD
    if (!found) blocked.include(visited);
#endif
    return found;
}

bool dfsPrune(Graph *g, int u, int start, HashSetStack& blocked, Path*& path, int tstart) {

    BlockedSet visited(g->getVertexNo());

    bool blck_flag = false;
    unsigned long vertex_visits = 0;
    bool found = dfsUtil(g, u, start, 0, blocked, visited, blck_flag, path, vertex_visits, tstart);

    vertexVisits.increment(vertex_visits);

#if defined(BLK_FORWARD) && defined(SUCCESSFUL_DFS_BLK)
    if (u != start && blck_flag) blocked.insert(u);
#endif

    if (!found) blocked.include(visited);

    return found;
}

/// The main backtracking function

namespace {

void followPath(Graph* g, EdgeData e, Cycle* current, HashSetStack& blocked,
                    CycleHist& result, Path* current_path, int tstart = -1);

void cyclesReadTarjan(Graph *g, EdgeData e, Cycle *current, HashSetStack &blocked,
                      CycleHist &result, pair<Path *, Path *> paths,
                      int tstart = -1) {

#ifdef NEIGHBOR_ORDER
    vector<Path*> allPaths;
#endif

    //exists a vertex w with a path to s*/
    for (int ind = g->offsArray[e.vertex]; ind < g->offsArray[e.vertex + 1]; ind++) {
        int w = g->edgeArray[ind].vertex;

        auto &tset = g->edgeArray[ind].tstamps;

        if (tstart != -1) {
            if (!edgeInTimeInterval(tstart, timeWindow, current->front(), e.vertex, tset)) continue;
        }

        if ((tstart == -1) && (w < current->front())) continue;
        if (blocked.exists(w)) continue;

        Path *current_path = NULL;
        if (paths.first && w == paths.first->back()) {
            current_path = paths.first; paths.first = NULL;
        } else if (paths.second && w == paths.second->back()) {
            current_path = paths.second; paths.second = NULL;
        } else {
            bool found = false;
            found = findPath(g, w, current->front(), blocked, current_path, tstart);
            if (!found) continue;
        }

#ifdef NEIGHBOR_ORDER
        allPaths.push_back(current_path);
    }

    for(auto current_path : allPaths) {
#endif

#ifdef BLK_FORWARD
//        HashSetStack new_blocked(blocked);
        blocked.incrementLevel();
        followPath(g, e, current, blocked, result, current_path, tstart);
        blocked.decrementLevel();
#else
        HashSetStack new_blocked(g->getVertexNo());
        for (auto c : *current) if (c != current->front()) new_blocked.insert(c);
        followPath(g, e, current, new_blocked, result, current_path, tstart);
#endif

    }
}

void followPath(Graph* g, EdgeData e, Cycle* current, HashSetStack& blocked, CycleHist& result, Path* current_path, int tstart) {

    Path *another_path = NULL;
    bool branching = false;

    int prev_vertex = -1;
    // Add vertices from the found path to the current path until a branching is found
    while(current_path->back() != current->front() && !branching) {
        prev_vertex = current_path->back();
        current_path->pop_back();
        current->push_back(prev_vertex);
        blocked.insert(prev_vertex);

        for(int ind = g->offsArray[prev_vertex]; ind < g->offsArray[prev_vertex+1]; ind++) {
            int u = g->edgeArray[ind].vertex;
            auto& tset = g->edgeArray[ind].tstamps;

            if(tstart != -1) {
                if(!edgeInTimeInterval(tstart, timeWindow, current->front(), prev_vertex, tset)) continue;
            }

            if(u != current_path->back() && ((tstart != -1) || (tstart == -1) && (u > current->front())) &&
              !blocked.exists(u) && !branching) {

                branching = dfsPrune(g, u, current->front(), blocked, another_path, tstart);

                if(branching) {
#ifndef BLK_FORWARD
                    delete another_path; another_path = NULL;
#endif
                    break;
                }
            }
        }
    }

    if(branching) {
#ifndef BLK_FORWARD
        if(current_path) delete current_path; current_path = NULL;
#endif
        cyclesReadTarjan(g, EdgeData(prev_vertex, -1), current, blocked, result, make_pair(current_path, another_path), tstart);
    }
    else {
        delete current_path; current_path = NULL;
        recordCycle(current, result);
    }

    // Remove all the vertices added after e.vertex
    while (1) {
        int back = current->back();
        if(back == e.vertex) break;
        current->pop_back();
    }
}

}

/// ************ coarse-grained Read-Tarjan algorithm with time window - top level ************

void allCyclesReadTarjanCoarseGrainedTW(Graph *g, CycleHist& result, int numThreads) {
    int process_rank = 0;
    int size_of_cluster = 1;

#ifdef MPI_IMPL
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif

    task_scheduler_init init(numThreads);
    tbb::combinable<CycleHist> pt_chist_jh;

    parallel_for(size_t(0), size_t(g->getVertexNo()), [&](size_t i) {
        if ((g->numNeighbors(i) != 0) && (g->numInEdges(i) != 0)) {
            parallel_for(size_t(g->offsArray[i]), size_t(g->offsArray[i + 1]), [&](size_t ind) {
                int w = g->edgeArray[ind].vertex;
                auto &tset = g->edgeArray[ind].tstamps;

                parallel_for(size_t(0), size_t(tset.size()), [&](size_t j) {

                    if((ind+j) % size_of_cluster == process_rank) {
                        auto &my_hist = pt_chist_jh.local();
                        int ts = tset[j];

                        HashSetStack blocked(g->getVertexNo());
                        Cycle *cycle = new Cycle();
                        cycle->push_back(i);

                        Path *current_path = NULL;
                        bool found = findPath(g, w, cycle->front(), blocked, current_path, ts);
                        if (found) followPath(g, EdgeData(i, -1), cycle, blocked, my_hist, current_path, ts);

                        delete cycle;
                    }
                });
            });
        }
    });

    combineCycleHistogram(pt_chist_jh, result);
}
