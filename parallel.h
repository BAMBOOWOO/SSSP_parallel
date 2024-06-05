#ifndef UTILS_PARALLEL_H
#define UTILS_PARALLEL_H

#include <iostream>
#include <vector>
using namespace std;

struct Edge {
    int to;
    int weight;
    
    Edge(int t, int w) : to(t), weight(w) {}
};

class Graph {
  
public:
    int V;
    vector<vector<pair<int,int>>> adjList;
    Graph(int vertices) : V(vertices) {
        adjList.resize(V);
    }

    void addEdge(int from, int to, int weight) {
        adjList[from].push_back({to, weight});
        adjList[to].push_back({from, weight}); // Assuming an undirected graph
    }

    void printGraph() {
        for (int v = 0; v < V; ++v) {
            cout << "Adjacency list of vertex " << v << "\n";
            for (auto edge: adjList[v]) {
                cout << "-> " << edge.first << " (weight: " << edge.second << ")\n";
            }
            cout << "\n";
        }
    }
};

void verify(const vector<int>& serial, const vector<int>& parallel, int V);
void printArr(const vector<int>& dist, int n, int size);
void printpath(const vector<int>& path, int size);
void parprintpath(const vector<int>& path, int size);
void BellmanFord(Graph graph, int src, vector<int>& serial, vector<int>& spath);
void BF_openMP(Graph graph, int src, vector<int>& parallel, vector<int>& ppath);
void BF_MPforMPI (const vector<vector<pair<int, int>>>& adjList, vector<int>& dist, vector<int>& ppath, int baserow, int gsize, int& updated, int rank, int num_procs);
void BF_MPI (Graph graph, int x, int y, int z, int gsize, vector<int>& parallel, vector<int>& ppath);
void flattenvec(const vector<vector<pair<int, int>>>& vec2d, vector<pair<int, int>>& vec1d);
void unflattenvec(vector<vector<pair<int, int>>>& vec2d, vector<pair<int, int>>& vec1d);
void preparedata(const vector<int>& localarray, vector<int>& sendarray, int rank, int flocalsize, int gsize);
#endif