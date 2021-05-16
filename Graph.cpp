#include "Graph.h"

using namespace std;

/*
 * ================================================================================================
 * Class Vertex Methods
 * ================================================================================================
 */

template <class T>
Vertex<T>::Vertex(T in): info(in) {}

template <class T>
void Vertex<T>::addEdge(Edge<T> *e) {
    adj.push_back(e);
    outgoing.push_back(e);
    e->dest->incoming.push_back(e);
}

template <class T>
bool Vertex<T>::operator<(Vertex<T> & vertex) const {
    return this->dist < vertex.dist;
}

template <class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template <class T>
vector<Edge<T> *>  Vertex<T>::getIncoming() const {
    return this->incoming;
}

template <class T>
vector<Edge<T> *>  Vertex<T>::getOutgoing() const {
    return this->outgoing;
}


/*
 * ================================================================================================
 * Class Edge Methods
 * ================================================================================================
 */

template <class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double cost, vector<T> complex_path): orig(o), dest(d), cost(cost){
    this->complex_path = complex_path;
}

template <class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double capacity, double cost, double flow):
        orig(o), dest(d), capacity(capacity), cost(cost), flow(flow){}

template <class T>
double  Edge<T>::getFlow() const {
    return this->flow;
}

template <class T>
vector<Edge<T> *>  Edge<T>::getComplexPath(){
    return this->complex_path;
}

template <class T>
void  Edge<T>::setComplexPath(vector<Edge<T> *> complex_path){
    this->complex_path = complex_path;
}


/*
 * ================================================================================================
 * Class Graph Methods
 * ================================================================================================
 */

template <class T>
Vertex<T> * Graph<T>::addVertex(const T &in) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in);
    vertexSet.push_back(v);
    return v;
}

template <class T>
Edge<T> * Graph<T>::addEdge(const T &sourc, const T &dest, double capacity, double cost, double flow) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d, capacity, cost, flow);
    s->addEdge(e);
    return e;
}

template <class T>
Edge<T> * Graph<T>::addEdge(const T &sourc, const T &dest, double cost, vector<T> complex_path) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d, cost, complex_path);
    s->addEdge(e);
    return e;
}

template <class T>
Vertex<T>* Graph<T>::findVertex(const T & inf) const {
    for (auto v : vertexSet)
        if (v->info == inf)
            return v;
    return nullptr;
}

template <class T>
double Graph<T>::getFlow(const T &sourc, const T &dest) const {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return 0.0;
    for (auto e : s->outgoing)
        if (e->dest == d)
            return e->flow;
    return 0.0;
}

template <class T>
vector<Vertex<T> *> Graph<T>::getVertexSet() const {
    return vertexSet;
}

/*
 * *************** DFS  ************
 */

template <class T>
std::vector<T> Graph<T>::dfs() const {
    std::vector<T> res;
    for(auto v : vertexSet){
        v->visited = false;
    }
    for(auto v : vertexSet){
        if(!v->visited)
            this->dfsVisit(v, res);
    }

    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent not yet visited, recursively.
 * Updates a parameter with the list of visited node contents.
 */
template <class T>
void Graph<T>::dfsVisit(Vertex<T> *v, std::vector<T> & res) const {
    v->visited = true;
    res.push_back(v->info);
    for(auto a: v->adj){
        if(!a.dest->visited)
            dfsVisit(a.dest, res);
    }
}

/*
 * *************** Shortest Path Problem  ***********
 */

template<class T>
void Graph<T>::dijkstraShortestPath(const T &origin) {
    for (Vertex<T>* vertex : vertexSet) {
        vertex->dist = INF;
        vertex->path = NULL;
        vertex->queueIndex = 0;
    }

    MutablePriorityQueue<Vertex<T>> vertexQueue;
    findVertex(origin)->dist = 0;
    vertexQueue.insert(findVertex(origin));

    while(!vertexQueue.empty()) {
        Vertex<T>* v = vertexQueue.extractMin();
        for (Edge<T>* edge : v->adj) {
            double oldDist = edge->getDest()->dist;
            if (edge->getDest()->dist > v->dist + edge->getCost()) {
                edge->getDest()->dist = v->dist + edge->getCost();
                edge->getDest()->path = v;
                if (oldDist == INF)
                    vertexQueue.insert(edge->dest);
                else
                    vertexQueue.decreaseKey(edge->dest);
            }
        }
    }
}

template<class T>
std::vector<T> Graph<T>::getPath(const T &origin, const T &dest) const{
    std::vector<T> res;

    Vertex<T>* v = this->findVertex(dest);
    if (v == nullptr || v->dist == INF) return res;

    while(true) {
        res.push_back(v->info);
        if (v->info == origin) break;
        v = v->path;
    }

    std::reverse(res.begin(), res.end());

    return res;
}

template <class T>
Graph<T> Graph<T>::generateInterestPointsGraph(vector<T> important_points){
    Graph<T> result;
    for(int i = 0; i < important_points.size(); i++){
        result.addVertex(important_points[i]);
    }
    for(int i = 0; i < important_points.size(); i++){
        T current_info = important_points[i];
        //Graph<T> copy = complete_graph;
        this->dijkstraShortestPath(current_info);
        //copy.dijkstraShortestPath(current_info);
        cout << endl << "I(" << i << ") -> currentInfo = " << current_info << endl;
        for(int j = 0; j < important_points.size(); j++){
            if(important_points[j] != current_info){
                vector<T> path = this->getPath(current_info, important_points[j]);
                double dist = this->findVertex(j)->getDist();
                cout << "J(" << j << ") -> Info = " << important_points[j] << endl;
                cout << "J(" << j << ") -> dist = " << dist << endl;
                cout << "J(" << j << ") -> path = ";
                for(int a = 0; a < path.size(); a++) cout << path[a] << ", ";
                cout << endl;
                result.addEdge(current_info, important_points[j], dist, path);
            }
        }
    }
    return result;
}


template<class T>
vector<Vertex<T> *> Graph<T>::nearestNeighbour(const T &origin_info) { //TODO: Test!
    Vertex<T> *origin = findVertex(origin_info);
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
