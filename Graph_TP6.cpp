#include "Graph_TP6.h"

using namespace std;

template<class T>
vector<Vertex<T> *> Graph<T>::nearestNeighbour(Vertex<T> *origin) { //TODO: Test!
    for (auto v : vertexSet)
        v->visited = false;
    vector<Vertex<T> *> path;
    origin->visited = true;
    path.push_back(origin);
    double minWeight;
    Vertex<T> *nextVertex = origin;
    do {
        minWeight = UINT_MAX; //TODO: get the maximum for a double
        auto adj = nextVertex->adj;
        for (auto edge : adj) {
            if (!edge.dest->visited && edge.weight < minWeight) {
                minWeight = edge.weight;
                nextVertex = edge.dest;
            }
        }
        nextVertex->visited = true;
        path.push_back(nextVertex);
    } while (minWeight != UINT_MAX);
    path.push_back(origin);
    return path;
}
