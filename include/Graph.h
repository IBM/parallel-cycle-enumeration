#ifndef GRAPH_H
#define GRAPH_H

#include <set>
#include <unordered_map>
#include "utils.h"
#include "Macros.h"

using namespace std;

struct EdgeData {
    int vertex; int tstamp; EdgeData(int v = 0, int t = 0) : vertex(v), tstamp(t) {}
};
typedef vector<int> TimestampSet;

struct ColElem {
    int vertex; TimestampSet tstamps; ColElem(int v = -1) : vertex(v) {}
};


class Graph {
public:
	Graph(){};
	~Graph() {delete [] offsArray; delete [] edgeArray; delete [] inOffsArray; delete [] inEdgeArray;};

    void readTemporalGraph(string path);
	int getVertexNo() { return vertexNo; }
    int getEdgeNo() { return edgeNo; }
    int numNeighbors(int node)  { return offsArray[node+1] - offsArray[node]; }
    int numInEdges(int node) { return inOffsArray[node+1] - inOffsArray[node]; }

	long vertexNo = 0, edgeNo = 0, tstampNo = 0;

	// Adjacency data, CSR format
	int* offsArray = nullptr;
    ColElem* edgeArray = nullptr;
    // Adjacency data, CSC format
    int* inOffsArray = nullptr;
    ColElem* inEdgeArray = nullptr;

};

inline void Graph::readTemporalGraph(string path) {
    ifstream graphFile(path);

    int edgeId = 0;
    unordered_map<int, map<int, TimestampSet>> adjacencyList;
    unordered_map<int, map<int, TimestampSet>> inEdgeList;

    int minTs = 2147483647, maxTs = 0;
    int maxNode = 0; tstampNo = 0;

    while(true) {
        string line; getline(graphFile, line);
        if(graphFile.eof()) break;
        if(line[0] == '%' || line[0] == '#') continue;

        stringstream ss(line);
        int first, second, third = 0;
        ss >> first >> second >> third;
        if(first != second) {
            tstampNo++;
            adjacencyList[first][second].push_back(third);
            inEdgeList[second][first].push_back(third);
            maxNode = max(maxNode, max(first, second));
            if(third != 0) minTs = min(minTs, third);
            maxTs = max(maxTs, third);
        }
        edgeId++;
    }
    vertexNo = maxNode + 1;

    for(auto pair : adjacencyList) edgeNo += pair.second.size();

    offsArray = new int [vertexNo + 2];
    edgeArray = new ColElem [edgeNo + 1];
    int currentOffs = 0; offsArray[0] = 0;
    for(int i = 0; i < vertexNo; i++) {
        if(adjacencyList.find(i) != adjacencyList.end()) {
            for(auto& pair : adjacencyList[i]) {
                sort(pair.second.begin(), pair.second.end());
                edgeArray[currentOffs].vertex = pair.first;
                edgeArray[currentOffs].tstamps = pair.second;
                currentOffs++;
            }
        }
        offsArray[i+1] = currentOffs;
    }

    inOffsArray = new int [vertexNo + 1];
    inEdgeArray = new ColElem [edgeNo];
    currentOffs = 0; inOffsArray[0] = 0;
    for(int i = 0; i < vertexNo; i++) {
        if(inEdgeList.find(i) != inEdgeList.end()) {
            for(auto& pair : inEdgeList[i]) {
                sort(pair.second.begin(), pair.second.end());
                inEdgeArray[currentOffs].vertex = pair.first;
                inEdgeArray[currentOffs].tstamps = pair.second;
                currentOffs++;
            }
        }
        inOffsArray[i+1] = currentOffs;
    }
    graphFile.close();
}
#endif