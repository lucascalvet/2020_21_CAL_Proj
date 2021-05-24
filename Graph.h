#ifndef PROJ_GRAPH_H_
#define PROJ_GRAPH_H_

#include <vector>
#include <queue>
#include <numeric>
#include <limits>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <graphviewer.h>
#include <sstream>
#include <string>
#include <cmath>
#include "MutablePriorityQueue.h"
#include "Utils.h"


constexpr auto INF = std::numeric_limits<double>::max();
constexpr auto UINF = std::numeric_limits<unsigned>::max();
constexpr double LATLNG_FACTOR = 50000.0;

class Van;

template<class T>
class Cluster;

template<class T>
class Vertex;

template<class T>
class Edge;

template<class T>
class Graph;

/*
 * ================================================================================================
 * Class Vertex
 * ================================================================================================
 */

template<class T>
class Vertex {
private:
    Vertex(T in, double lat, double lng, unsigned hour, unsigned tolerance, unsigned max_tolerance,
           unsigned quantity = 0);

    Vertex(T in, double lat, double lng, unsigned quantity = 0);

    /// unique identifier of the vertex
    T info;
    /// latitude (or x) position of the vertex
    double lat;
    /// longitude (or y) position of the vertex
    double lng;
    /// outgoing edges of the vertex
    std::vector<Edge<T> *> outgoing;
    /// incoming edges of the vertex
    std::vector<Edge<T> *> incoming;
    /// all edges of the vertex
    std::vector<Edge<T> *> adj;
    /// used for path finding
    bool visited = false;
    Vertex<T> *path = NULL;
    double dist = INF;
    /// required by MutablePriorityQueue
    int queueIndex = 0;
    /// hour that the client wants the delivery (minutes from midnight)
    unsigned hour = 0;
    /// tolerance for the delivery, after the hour, from which a delay is counted
    unsigned tolerance = 0;
    /// maximum tolerance for the delivery
    unsigned max_tolerance = 0;
    /// quantity of breads in the order
    unsigned quantity = 0;
    /// used for path finding, the hour that the client was visited
    unsigned visited_at;

    void addEdge(Edge<T> *e);

    bool operator<(Vertex<T> &vertex) const; // required by MutablePriorityQueue

public:
    T getInfo() const;

    double getDist() { return dist; }

    std::vector<Edge<T> *> getIncoming() const;

    std::vector<Edge<T> *> getOutgoing() const;

    double getLat() const { return lat; }

    void setLat(double lat) { Vertex::lat = lat; }

    double getLng() const { return lng; }

    void setLng(double lng) {
        Vertex::lng = lng;
    }

    unsigned int getHour() const {
        return hour;
    }

    void setHour(unsigned int hour) {
        Vertex::hour = hour;
    }

    unsigned int getTolerance() const {
        return tolerance;
    }

    void setTolerance(unsigned int tolerance) {
        Vertex::tolerance = tolerance;
    }

    unsigned int getMaxTolerance() const {
        return max_tolerance;
    }

    void setMaxTolerance(unsigned int maxTolerance) {
        max_tolerance = maxTolerance;
    }

    unsigned int getQuantity() const {
        return quantity;
    }

    void setQuantity(unsigned int quantity) {
        Vertex::quantity = quantity;
    }

    unsigned int getVisitedAt() const {
        return visited_at;
    }

    void setVisitedAt(unsigned int visitedAt) {
        visited_at = visitedAt;
    }

    void setTimes(unsigned hour, unsigned tolerance, unsigned max_tolerance) {
        this->hour = hour;
        this->tolerance = tolerance;
        this->max_tolerance = max_tolerance;
    }

    void removeEdge(Edge<T> *e);

    friend class Graph<T>;

    friend class MutablePriorityQueue<Vertex<T>>;
};

template<class T>
Vertex<T>::Vertex(T in, double lat, double lng, unsigned quantity): info(in), lat(lat), lng(lng), quantity(quantity) {}

template<class T>
Vertex<T>::Vertex(T in, double lat, double lng, unsigned hour, unsigned tolerance, unsigned max_tolerance,
                  unsigned quantity): info(in), lat(lat), lng(lng), hour(hour), tolerance(tolerance),
                                      max_tolerance(max_tolerance), quantity(quantity) {}

/**
 * Add outgoing edge to a vertex
 *
 * @param e pointer to the new edge
 */
template<class T>
void Vertex<T>::addEdge(Edge<T> *e) {
    adj.push_back(e);
    outgoing.push_back(e);
    e->dest->incoming.push_back(e);
    e->dest->adj.push_back(e);
}

/**
 * Remove outgoing edge from a vertex
 *
 * @param e pointer to the edge to be removed
 */
template<class T>
void Vertex<T>::removeEdge(Edge<T> *e) {
    for (auto edge = adj.begin(); edge < adj.end(); edge++) {
        if (*edge == e) {
            adj.erase(edge);
            break;
        }
    }
    for (auto edge = outgoing.begin(); edge < outgoing.end(); edge++) {
        if (*edge == e) {
            outgoing.erase(edge);
            break;
        }
    }
    for (auto edge = e->dest->incoming.begin(); edge < e->dest->incoming.end(); edge++) {
        if (*edge == e) {
            e->dest->incoming.erase(edge);
            break;
        }
    }
    delete e;
}

/**
 * Compare two vertices (required by MutablePriorityQueue)
 *
 * @param vertex the vertex to compare to
 * @return if the vertex dist is lower
 */
template<class T>
bool Vertex<T>::operator<(Vertex<T> &vertex) const {
    return this->dist < vertex.dist;
}

/**
 * Get the info of a Vertex
 *
 * @return the info of a Vertex
 */
template<class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

/**
 * Get the incoming edges of a Vertex
 *
 * @return a vector of pointers to the incoming edges of a Vertex
 */
template<class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

/**
 * Get the outgoing edges of a Vertex
 *
 * @return a vector of pointers to the outgoing edges of a Vertex
 */
template<class T>
std::vector<Edge<T> *> Vertex<T>::getOutgoing() const {
    return this->outgoing;
}

/*
 * ================================================================================================
 * Class Edge
 * ================================================================================================
 */

template<class T>
class Edge {
    /// the origin vertex
    Vertex<T> *orig;
    /// the destination vertex
    Vertex<T> *dest;
    /// the underlying path of an edge (for interest points graphs)
    std::vector<T> complex_path;
    /// the edge cost
    double cost;
    /// required to view Van Paths in GraphViewer
    int vanID = 0;

public:
    Edge(Vertex<T> *o, Vertex<T> *d, double cost, std::vector<T> complex_path);

    Edge(Vertex<T> *o, Vertex<T> *d, double cost);

    Edge(Vertex<T> *o, Vertex<T> *d, bool lat_long);

    Vertex<T> *getOrigin() { return orig; }

    Vertex<T> *getDest() { return dest; }

    int getVanID() { return vanID; }

    void setVanID(int vanID) { this->vanID = vanID; }

    double getCost() const { return cost; }

    void setCost(double cost) { this->cost = cost; }

    //std::vector<Edge<T> *> getComplexPath();
    std::vector<T> getComplexPath() {
        return this->complex_path;
    }

    void setComplexPath(std::vector<Edge<T> *> complex_path) {
        this->complex_path = complex_path;
    }

    friend class Graph<T>;

    friend class Vertex<T>;
};

template<class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double cost, std::vector<T> complex_path): orig(o), dest(d), cost(cost) {
    this->complex_path = complex_path;
}

template<class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double cost):
        orig(o), dest(d), cost(cost) {}

template<class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, bool lat_long):
        orig(o), dest(d) {
    if (lat_long)
        cost = calculateDist(this->orig->getLng(), this->orig->getLat(), this->dest->getLng(), this->dest->getLat());
    else cost = calculateDist(this->orig->getLat(), this->orig->getLng(), this->dest->getLat(), this->dest->getLng());
}

template<class T>
class ClientCompare {
public:
    bool operator()(const Vertex<T> *v1, const Vertex<T> *v2) {
        return v1->getHour() + v1->getMaxTolerance() > v2->getHour() + v2->getMaxTolerance();;
    }
};

template<class T>
class ClientArrivalCompare {
public:
    bool operator()(Vertex<T> *v1, Vertex<T> *v2) {
        return v1->getHour() + v1->getTolerance() > v2->getHour() + v2->getTolerance();;
    }
};

/*
 * ================================================================================================
 * Class Graph
 * ================================================================================================
 */

template<class T>
class Graph {
private:
    /// list of vertices in the graph
    std::vector<Vertex<T> *> vertexSet;
    /// indicates if it is an interest points graph, suitable for the corresponding algorithms
    bool ipGraph = false;
    bool latLng = true;

    unsigned start_time = 0;
    double weight = 0.5;
    unsigned visit_time = 5;
    unsigned early_time = 10;
    unsigned velocity = 50;

public:
    //~Graph();

    void nearestNeighbour(const T &origin_info);

    void nearestNeighbourTimes(const T &origin_info);

    void dijkstraShortestPath(const T &origin);

    void heldKarp(const T &origin);

    double bruteForceTimes(const T &origin);

    std::pair<double, std::vector<Vertex<T> *>>
    bruteForceTimesRec(Vertex<T> *origin,
                       std::priority_queue<Vertex<T> *, std::vector<Vertex<T> *>, ClientCompare<T>>

                       client_queue,
                       std::vector<Vertex<T> *> path = std::vector<Vertex<T> *>()
    );

    std::vector<std::vector<T>> dfsConnectivity() const;

    bool checkConectivity(std::vector<T> ids);

    std::vector<T> dfs() const;

    void dfsVisit(Vertex<T> *v, std::vector<T> &res, bool reverse = false) const;

    Graph<T> generateInterestPointsGraph(std::vector<T> important_points);

    Vertex<T> *findVertex(const T &inf) const;

    Edge<T> *findEdge(const T &og, const T &dest) const;

    Edge<T> *findEdge(const Vertex<T> *og, const Vertex<T> *dest) const;

    std::vector<Vertex<T> *> getVertexSet() const;

    unsigned int getStartTime() const;

    void setStartTime(unsigned int startTime);

    unsigned int getVisitTime() const;

    bool getLatLng() const { return this->lat_lng; }

    void setLatLng(bool lat_lng) { this->lat_lng = lat_lng; }

    void setVisitTime(unsigned int visitTime);

    unsigned int getEarlyTime() const;

    void setEarlyTime(unsigned int earlyTime);

    unsigned getTotalQuantity() const {unsigned totalQuantity = 0; for(Vertex<T> * v : vertexSet) totalQuantity += v->getQuantity(); return totalQuantity;};

    Vertex<T> *addVertex(const T &in);

    //Vertex<T> *addVertex(const T &in, const double &x, const double &y);

    Vertex<T> *addVertex(const T &in, const double &lat, const double &lng, const unsigned &quantity = 0);

    Vertex<T> *
    addVertex(const T &in, const double &lat, const double &lng, const unsigned &hour, const unsigned &tolerance,
              const unsigned &max_tolerance, const unsigned &quantity = 0);

    //Vertex<T> *addVertex(const T &in, const double &x, const double &y, const unsigned &hour, const unsigned &tolerance, const unsigned &max_tolerance);

    Edge<T> *addEdge(const T &sourc, const T &dest, double capacity, double cost, double flow = 0);

    Edge<T> *addEdge(const T &sourc, const T &dest, double cost, std::vector<T> complex_path);

    Edge<T> *addEdge(const T &sourc, const T &dest, double cost);

    Edge<T> *addEdge(const T &sourc, const T &dest, bool lat_lng);

    std::pair<std::vector<T>, double> getPath(const T &origin, const T &dest) const;

    unsigned getVelocity() { return this->velocity; }

    void setVelocity(unsigned velocity) { this->velocity = velocity; }

    void viewGraph();

    void viewGraphIP(const Graph<T> &igraph);

    void viewGraphPath(std::vector<T> path, std::vector<T> interest_points, bool remove_extra_nodes = false,
                       bool remove_extra_edges = false, bool show_w = false, bool show_nid = false, bool multiple_vans = true);

    void viewGraphPathIP(Graph<T> igraph, std::vector<T> path, bool remove_extra_nodes = false,
                         bool remove_extra_edges = false, bool show_w = false, bool show_nid = false, bool multiple_vans = true);

    void importGraph(std::string vertex_filename, std::string edges_filename, bool lat_lng);

    void printTimes();

    std::vector<T> getOverlapClients(const T &info);

    std::vector<T> getOverlapClientsTravelling(const T &info);

    std::vector<Vertex<T> *> getOverlapClientsTravelling(Vertex<T> *info);

    std::vector<T> getPerfectOverlapClientsTravelling(const T &info);

    std::vector<T> getPerfectOverlapClientsTravelling(std::vector<T> remaining_clients, const T &info);

    double costFunctionStep(T og, std::vector<T> path, T new_element);

    double costFunctionTotal(T og, std::vector<T> path);

    double costFunctionTotal(Vertex<T> *og, std::vector<Vertex<T> *> path);

    std::vector<Cluster<T>> getClusters(const T &origin);

    std::vector<std::pair<Van, std::vector<Vertex<T> *>>> dividingClustersBrute(std::vector<Van> vans, const T &origin);

    std::vector<std::pair<Van, std::vector<Vertex<T> *>>>
    dividingClustersGreedy(std::vector<Van> vans, const T &origin);

    void followVansPath(std::vector<std::pair<Van, std::vector<Vertex<T> *>>> paths_dist, const T &og);

    std::vector<T> vanPathtoViewPath(std::vector<std::pair<Van, std::vector<Vertex<T > *>>> van_path, const T &og);

};

/*
template<class T>
Graph<T>::~Graph() {
    for (auto v : vertexSet) {
        for (auto out : v->outgoing) {
            delete out;
        }
        delete v;
    }
}
 */

template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in);
    vertexSet.push_back(v);
    return v;
}

/*
template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in, const double &x, const double &y) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in, x, y);
    vertexSet.push_back(v);
    return v;
}
 */

template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in, const double &lat, const double &lng, const unsigned &quantity) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in, lat, lng, quantity);
    vertexSet.push_back(v);
    return v;
}

/*
template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in, const double &x, const double &y, const unsigned &hour, const unsigned &tolerance, const unsigned &max_tolerance){
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in, x, y, hour, tolerance, max_tolerance);
    vertexSet.push_back(v);
    return v;
}
 */

template<class T>
Vertex<T> *
Graph<T>::addVertex(const T &in, const double &lat, const double &lng, const unsigned &hour, const unsigned &tolerance,
                    const unsigned &max_tolerance, const unsigned &quantity) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in, lat, lng, hour, tolerance, max_tolerance, quantity);
    vertexSet.push_back(v);
    return v;
}


template<class T>
Edge<T> *Graph<T>::addEdge(const T &sourc, const T &dest, double capacity, double cost, double flow) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d, capacity, cost, flow);
    s->addEdge(e);
    return e;
}

template<class T>
Edge<T> *Graph<T>::addEdge(const T &sourc, const T &dest, double cost, std::vector<T> complex_path) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d, cost, complex_path);
    s->addEdge(e);
    return e;
}

template<class T>
Edge<T> *Graph<T>::addEdge(const T &sourc, const T &dest, double cost) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d, cost);
    s->addEdge(e);
    return e;
}

template<class T>
Edge<T> *Graph<T>::addEdge(const T &sourc, const T &dest, bool lat_lng) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d, lat_lng);
    s->addEdge(e);
    return e;
}

template<class T>
Vertex<T> *Graph<T>::findVertex(const T &inf) const {
    for (auto v : vertexSet)
        if (v->info == inf)
            return v;
    return nullptr;
}

template<class T>
Edge<T> *Graph<T>::findEdge(const T &og, const T &dest) const {
    for (auto v : vertexSet) {
        for (auto edge : v->getOutgoing()) {
            if (edge->orig->getInfo() == og && edge->dest->getInfo() == dest) return edge;
        }
    }
    return nullptr;
}

template<class T>
Edge<T> *Graph<T>::findEdge(const Vertex<T> *og, const Vertex<T> *dest) const {
    for (auto out : og->outgoing) {
        if (out->dest == dest)
            return out;
    }
    return nullptr;
}

template<class T>
std::vector<Vertex<T> *> Graph<T>::getVertexSet() const {
    return vertexSet;
}

template<class T>
unsigned int Graph<T>::getStartTime() const {
    return start_time;
}

template<class T>
void Graph<T>::setStartTime(unsigned int startTime) {
    start_time = startTime;
}

template<class T>
unsigned int Graph<T>::getVisitTime() const {
    return visit_time;
}

template<class T>
void Graph<T>::setVisitTime(unsigned int visitTime) {
    visit_time = visitTime;
}

template<class T>
unsigned int Graph<T>::getEarlyTime() const {
    return early_time;
}

template<class T>
void Graph<T>::setEarlyTime(unsigned int earlyTime) {
    early_time = earlyTime;
}

/*
 * ================================================================================================
 * DFS
 * ================================================================================================
 */

template<class T>
std::vector<std::vector<T>> Graph<T>::dfsConnectivity() const {
    std::vector<T> visit_order = dfs();
    for (auto v : vertexSet) {
        v->visited = false;
    }
    Vertex<T> *v;
    std::vector<T> scc;
    std::vector<std::vector<T>> res;
    for (int i = visit_order.size() - 1; i >= 0; i--) {
        v = findVertex(visit_order.at(i));
        if (!v->visited) {
            dfsVisit(v, scc, true);
            res.push_back(scc);
            scc.clear();
        }
    }
    return res;
}

template<class T>
std::vector<T> Graph<T>::dfs() const {
    std::vector<T> res;
    for (auto v : vertexSet) {
        v->visited = false;
    }
    for (auto v : vertexSet) {
        if (!v->visited)
            dfsVisit(v, res);
    }
    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent not yet visited, recursively.
 * Updates a parameter with the list of visited node contents.
 */
template<class T>
void Graph<T>::dfsVisit(Vertex<T> *v, std::vector<T> &res, bool reverse) const {
    v->visited = true;
    std::vector<Edge<T> *> edges;
    if (reverse) {
        for (auto e : v->incoming) {
            if (!e->orig->visited)
                dfsVisit(e->orig, res, reverse);
        }
    } else {
        for (auto e : v->outgoing) {
            if (!e->dest->visited)
                dfsVisit(e->dest, res, reverse);
        }
    }
    res.push_back(v->info);
}

/*
 * ================================================================================================
 * Shortest Path Problem
 * ================================================================================================
 */

template<class T>
void Graph<T>::dijkstraShortestPath(const T &origin) {
    for (Vertex<T> *vertex : vertexSet) {
        vertex->dist = INF;
        vertex->path = NULL;
        vertex->queueIndex = 0;
    }

    MutablePriorityQueue<Vertex<T>> vertexQueue;
    findVertex(origin)->dist = 0;
    vertexQueue.insert(findVertex(origin));

    while (!vertexQueue.empty()) {
        Vertex<T> *v = vertexQueue.extractMin();
        for (Edge<T> *edge : v->outgoing) {
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
std::pair<std::vector<T>, double> Graph<T>::getPath(const T &origin, const T &dest) const {
    std::vector<T> path;
    double total_dist = 0;

    Vertex<T> *v = this->findVertex(dest);
    if (v == nullptr || v->dist == INF) return std::pair<std::vector<T>, double>(path, total_dist);

    while (true) {
        path.push_back(v->info);
        total_dist += v->dist;
        v = v->path;
        if (v == nullptr) {
            return std::pair<std::vector<T>, double>(path, total_dist);
        }
        if (v->info == origin) {
            path.push_back(v->info);
            break;
        }
    }

    std::reverse(path.begin(), path.end());

    return std::pair<std::vector<T>, double>(path, total_dist);
}

template<class T>
Graph<T> Graph<T>::generateInterestPointsGraph(std::vector<T> important_points) {
    Graph<T> result;
    if (important_points.empty()) return result;
    result.visit_time = this->visit_time;
    result.early_time = this->early_time;
    result.start_time = this->start_time;
    result.velocity = this->velocity;
    std::cout << "Checking Con..." << std::endl;
    if (!checkConectivity(important_points)) {
        std::cout << "Conectivity Invalid!!!" << std::endl;
        return result;
    }
    std::cout << "Con Good!" << std::endl;

    for (int i = 0; i < important_points.size(); i++) {
        Vertex<T> *v = findVertex(important_points[i]);
        result.addVertex(important_points[i], v->getLat(), v->getLng(), v->getHour(), v->getTolerance(),
                         v->getMaxTolerance(), v->getQuantity());
    }

    if (important_points.size() == 1) return result;
    for (int i = 0; i < important_points.size(); i++) {
        T current_info = important_points[i];
        //Graph<T> copy = complete_graph;
        this->dijkstraShortestPath(current_info);
        //copy.dijkstraShortestPath(current_info);
        for (int j = 0; j < important_points.size(); j++) {
            if (important_points[j] != current_info) {
                std::vector<T> path = (getPath(current_info, important_points[j])).first;
                double dist = findVertex(important_points[j])->getDist();
                result.addEdge(current_info, important_points[j], dist, path);
            }
        }
    }
    return result;
}

class Van {
private:
    unsigned id;
    unsigned visit_time;
    unsigned capacity;

public:
    Van(unsigned id, unsigned visit_time, unsigned capacity=1) {
        this->id = id;
        this->visit_time = visit_time;
        this->capacity = capacity;
    }

    unsigned getID() const { return this->id; }

    unsigned getVisitTime() const { return this->visit_time; }

    void setVisitTime(unsigned visit_time) { this->visit_time = visit_time; }

    unsigned getCapacity() const { return this->capacity; }

    void setCapacity(unsigned capacity) { this->capacity = capacity; }

};

class VanCompare {
public:
    bool operator()(Van v1, Van v2) {
        return v1.getVisitTime() > v2.getVisitTime();
    }
};

template<class T>
class Cluster {
private:
    Graph<T> * graph;
    std::vector<T> info_set;
    std::vector<Vertex<T> *> vertex_set;

public:
    Cluster(Graph<T> * graph, std::vector<T> info_set) {
        this->graph = graph;
        this->info_set = info_set;
    }

    unsigned getSize() { return info_set.size(); }

    std::vector<T> getInfoSet() { return info_set; }

    std::vector<T> getVertexSet() { return vertex_set; }

    void addInfo(T info) { info_set.push_back(info); }

    void addVertex(Vertex<T> *v) {
        vertex_set.push_back(v);
        addInfo(v->info);
    }

    void PrintInfo(){
        for(T info : info_set){
            std::cout << info << " ";
        }
        std::cout << "END" << std::endl;
    }

    unsigned getTotalQuantity() const {
        unsigned totalQuantity = 0;
        for (T info : info_set){
            Vertex<T> * v = graph->findVertex(info);
            if(v == nullptr) return UINF;
            totalQuantity += v->getQuantity();
        }
        return totalQuantity;
    };

    bool isInSet(T new_info) {
        for (T info : info_set) {
            if (new_info == info) return true;
        }
        return false;
    }

    std::pair<Vertex<T> *, unsigned> getArrivalHour() {
        unsigned minHour = UINF;
        Vertex<T> *selected_vertex = nullptr;
        for (T info : info_set) {
            Vertex<T> *vertex = graph->findVertex(info);
            if (vertex != nullptr) {
                unsigned vertex_mh = vertex->getHour() + vertex->getTolerance();
                if (vertex_mh < minHour) {
                    minHour = vertex_mh;
                    selected_vertex = vertex;
                }
            }
        }
        return std::pair<Vertex<T> *, unsigned>(selected_vertex, minHour);
    };

    std::pair<Vertex<T> *, unsigned> getDepartureHour() {
        unsigned maxHour = 0;
        Vertex<T> *selected_vertex = nullptr;
        for (T info : info_set) {
            Vertex<T> *vertex = graph->findVertex(info);
            if (vertex != nullptr) {
                //unsigned vertex_mh = vertex->getHour() + vertex->getTolerance() + graph.getVisitTime();
                unsigned vertex_mh = vertex->getHour() + vertex->getMaxTolerance() + graph->getVisitTime();
                if (vertex_mh > maxHour) {
                    maxHour = vertex_mh;
                    selected_vertex = vertex;
                }
            }
        }
        return std::pair<Vertex<T> *, unsigned>(selected_vertex, maxHour);
    };

    Vertex<T> *getArrivalVertex() {
        return getArrivalHour().first;
    }

    Vertex<T> *getDepartureVertex() {
        return getDepartureHour().first;
    }

    unsigned getMaxArrivalHour() {
        Vertex<T> *v = getArrivalHour().first;
        return v->getHour() + v->getMaxTolerance();
    }

    double calculateClusterCost(Cluster<T> otherCluster);

    double calculateOriginCost(const T &origin, bool to_origin = false);

    void mergeCluster(Cluster<T> otherCluster);

    bool hasCommon(Cluster<T> otherCluster);

    std::vector<Vertex<T> *> exportVertices(){
        std::priority_queue<Vertex<T> *, std::vector<Vertex<T> *>, ClientArrivalCompare<T>> vertex_queue;
        std::vector<Vertex<T> *> result;
        for(T info: info_set){
            Vertex<T> * v = graph->findVertex(info);
            if(v == nullptr) return {};
            vertex_queue.push(v);
        }
        while(!vertex_queue.empty()){
            result.push_back(vertex_queue.top());
            vertex_queue.pop();
        }
        return result;
    }

    /*
    Cluster<T> operator=(Cluster<T> otherCluster){
        //this->info_set = otherCluster.info_set;
        //return this;
        return Cluster<T>(graph, otherCluster.getInfoSet());
    }
     */
};


template<class T>
class ClusterCompareArrival {
public:
    bool operator()(Cluster<T> c1, Cluster<T> c2) {
        return c1.getArrivalHour().second > c2.getArrivalHour().second;
    }
};

template<class T>
class ClusterCompareDeparture {
public:
    bool operator()(Cluster<T> c1, Cluster<T> c2) {
        return c1.getDepartureHour().second > c2.getDepartureHour().second;
    }
};

template<class T>
class ClusterCompareSize {
public:
    bool operator()(Cluster<T> c1, Cluster<T> c2) {
        return c1.getSize() < c2.getSize();
    }
};

template<class T>
void Cluster<T>::mergeCluster(Cluster<T> otherCluster) {
    for (T info : otherCluster.info_set) {
        if (!isInSet(info)) addInfo(info);
    }
}

template<class T>
bool Cluster<T>::hasCommon(Cluster<T> otherCluster) {
    for (T info : info_set) {
        for (T other_info : otherCluster.getInfoSet()) {
            if (other_info == info) return true;
        }
    }
    return false;
}

template<class T>
double Cluster<T>::calculateClusterCost(Cluster<T> otherCluster) {
    std::pair<Vertex<T> *, unsigned> departure = getDepartureHour();
    std::pair<Vertex<T> *, unsigned> arrival = otherCluster.getArrivalHour();
    if (arrival.first == nullptr || departure.first == nullptr) return INF;

    std::cout << "d->" << departure.second << " vs a->" << arrival.second << std::endl;

    Edge<T> *edge = graph->findEdge(departure.first, arrival.first); //TODO: removed getInfo()
    if (edge == nullptr) return INF;

    double travelling_time = edge->getCost() / graph->getVelocity();

    std::cout << "d+t->" << departure.second + travelling_time << " vs(<= 0) a->" << arrival.second << std::endl;
    std::cout << "d+t->" << departure.second + travelling_time << " vs(> INF) max_a->" << otherCluster.getMaxArrivalHour() << std::endl;
    std::cout << "return delay->" << departure.second + travelling_time - arrival.second << std::endl;

    if (departure.second + travelling_time <= arrival.second) return 0;
    else if (departure.second + travelling_time > otherCluster.getMaxArrivalHour()) return INF;
    else return departure.second + travelling_time - arrival.second;
}


template<class T>
double Cluster<T>::calculateOriginCost(const T &origin, bool to_origin) {
    if (to_origin) {
        Edge<T> *edge = graph->findEdge(getDepartureVertex()->getInfo(), origin);
        if (edge == nullptr) return INF;
        return edge->getCost();
    } else {
        std::pair<Vertex<T> *, unsigned> arrival = getArrivalHour();
        Edge<T> *edge = graph->findEdge(origin, arrival.first->getInfo());
        if (edge == nullptr) return -1;
        double delay = graph->getStartTime() + edge->getCost() - arrival.second;
        if (delay < 0) return 0;
        if (graph->getStartTime() + edge->getCost() > getMaxArrivalHour()) return INF;
        return delay;
    }
}

template<class T>
void Graph<T>::heldKarp(const T &origin) {
    //Initialize vertices vector with origin in last place
    Vertex<T> *orig = findVertex(origin);
    if (orig == nullptr) {
        std::cout << "[Held-Karp] Invalid origin!\n";
        return;
    }
    std::vector<Vertex<T> *> vertices = vertexSet;
    if (vertices.size() == 1) {
        std::cout << "[Held-Karp] Graph does not have clients!\n";
    }

    vertices.erase(std::find(vertices.begin(), vertices.end(), orig));
    vertices.push_back(orig);

    //Build distance matrix
    std::cout << "[Held-Karp] Building distance matrix...\n";
    std::vector<std::vector<double>> distance(vertices.size(), std::vector<double>(vertices.size(), 0));
    std::vector<Edge<T> *> edges;
    bool found_edge = false;
    for (unsigned i = 0; i < vertices.size(); i++) {
        edges = vertices.at(i)->outgoing;
        for (unsigned j = 0; j < vertices.size(); j++) {
            if (i == j) continue;
            for (auto edge = edges.begin(); edge < edges.end(); edge++) {
                if ((*edge)->dest == vertices.at(j)) {
                    distance.at(i).at(j) = (*edge)->cost;
                    edges.erase(edge);
                    found_edge = true;
                    break;
                }
            }
            if (!found_edge) {
                std::cout << "[Held-Karp] No edge connecting " << vertices.at(i)->info << " to "
                          << vertices.at(j)->info << "!\n";
                return;
            }
        }
    }

    std::cout << "[Held-Karp] Calculate shortest paths...\n";
    //Initialize best distance and corresponding path vectors
    std::vector<std::vector<unsigned>> best(1 << (vertices.size() - 1),
                                            std::vector<unsigned>(vertices.size(), UINT_MAX));
    std::vector<std::vector<unsigned>> path(1 << (vertices.size() - 1),
                                            std::vector<unsigned>(vertices.size(), UINT_MAX));
    unsigned temp_dist;
    //Iterate through combinations of vertices encoded in bits
    for (unsigned visited = 1; visited < (1 << (vertices.size() - 1)); visited++) {
        //Iterate through the possible last vertices in a path
        for (unsigned last = 0; last < vertices.size() - 1; last++) {
            //The last vertices for a certain visited combination need to be part of that combination
            if (!(visited & 1 << last)) continue;

            //If there is only one vertex in visited, the minimum distance is the one from the origin to last
            if (visited == 1 << last) {
                best.at(visited).at(last) = distance.at(vertices.size() - 1).at(last);
                path.at(visited).at(last) = vertices.size() - 1;
            } else {
                unsigned prev_visited = visited ^(1 << last); //possible previous visited indexes
                for (unsigned cand_prev = 0; cand_prev < vertices.size() - 1; cand_prev++) {
                    if (!(prev_visited & 1 << cand_prev)) continue; //discard impossible previous indexes
                    temp_dist = best.at(prev_visited).at(cand_prev) + distance.at(cand_prev).at(last);
                    if (temp_dist < best.at(visited).at(last)) {
                        best.at(visited).at(last) = temp_dist;
                        path.at(visited).at(last) = cand_prev;
                    }
                }
            }
        }
    }

    std::cout << "[Held-Karp] Calculate final shortest path...\n";
    //Get the cheapest path from the precomputed path lengths using all but the origin vertex
    unsigned final_dist = UINF, second_to_last = 0;
    for (unsigned last = 0; last < vertices.size() - 1; last++) {
        temp_dist = best.at((1 << (vertices.size() - 1)) - 1).at(last) + distance.at(last).at(vertices.size() - 1);
        if (temp_dist < final_dist) {
            final_dist = temp_dist;
            second_to_last = last;
        }
    }

    std::cout << "[Held-Karp] Populate path and dist attributes...\n";
    //Populate the path and dist attribute in the vertices
    unsigned last = second_to_last, visited = (1 << (vertices.size() - 1)) - 1, temp_last;
    vertices.at(vertices.size() - 1)->path = vertices.at(second_to_last);
    vertices.at(vertices.size() - 1)->dist = distance.at(second_to_last).at(vertices.size() - 1);
    for (unsigned i = 0; i < vertices.size() - 1; i++) {
        vertices.at(last)->path = vertices.at(path.at(visited).at(last));
        vertices.at(last)->dist = distance.at(path.at(visited).at(last)).at(last);
        temp_last = path.at(visited).at(last);
        visited ^= (1 << last);
        last = temp_last;
    }
    if (visited == 0) {
        std::cout << "[Held-Karp] All good!\n";
    } else std::cout << "[Held-Karp] All is not good :)!\n";
}

template<class T>
double Graph<T>::bruteForceTimes(const T &origin) {
    Vertex<T> *orig = findVertex(origin);
    if (orig == nullptr) {
        std::cout << "[Brute Force] Invalid origin!\n";
        return -1;
    }

    std::priority_queue<Vertex<T> *, std::vector<Vertex<T> *>, ClientCompare<T>> client_queue;

    for (Vertex<T> *v : vertexSet) {
        if (v != orig)
            client_queue.push(v);
    }

    std::vector<Vertex<T> *> path;

    std::pair<double, std::vector<Vertex<T> *>> res = bruteForceTimesRec(orig, client_queue, path);

    if (res.first == INF) return res.first;

    path = res.second;
    if (path.size() > 0)
        res.second.at(0)->path = orig;
    for (int i = 1; i < path.size(); i++) {
        path.at(i)->path = path.at(i - 1);
    }
    if (path.size() > 0)
    orig->path = path.at(path.size() - 1);

    unsigned current_time = start_time;
    Edge<T> *edge;

    if (path.size() > 0) {
        edge = findEdge(orig, path.at(0));
        if (edge == nullptr) {
            std::cout << "[Brute Force Times] Couldn't find edge!\n";
            return res.first;
        } else {
            current_time += edge->cost / velocity;
            if (current_time < edge->dest->hour - early_time) { //arrived early
                current_time = edge->dest->hour - early_time;
            }
            path.at(0)->visited_at = current_time;
            current_time += visit_time;
            path.at(0)->dist = edge->cost;
        }
    }

    for (int i = 1; i < path.size(); i++) {
        edge = findEdge(orig, res.second.at(0));
        current_time += edge->cost / velocity;
        if (current_time < edge->dest->hour - early_time) { //arrived early
            current_time = edge->dest->hour - early_time;
        }
        path.at(0)->visited_at = current_time;
        current_time += visit_time;
        path.at(0)->dist = edge->cost;
    }

    return res.first;
}

template<class T>
std::pair<double, std::vector<Vertex<T> *>>
Graph<T>::bruteForceTimesRec(Vertex<T> *origin,
                             std::priority_queue<Vertex<T> *, std::vector<Vertex<T> *>, ClientCompare<T>> client_queue,
                             std::vector<Vertex<T> *> path) {
    if (client_queue.empty())
        return std::pair<double, std::vector<Vertex<T> *>>(costFunctionTotal(origin, path), path);
    Vertex<T> *next_min = client_queue.top();
    std::vector<Vertex<T> *> overlapping = getOverlapClientsTravelling(next_min);
    std::vector<Vertex<T> *> temp_clients;
    std::vector<Vertex<T> *> temp_path;
    std::priority_queue<Vertex<T> *, std::vector<Vertex<T> *>, ClientCompare<T>> temp_queue;
    double min_cost = INF;
    std::vector<Vertex<T> *> min_path;
    std::pair<double, std::vector<Vertex<T> *>> res;

    temp_path = path;

    for (auto client : overlapping) {
        temp_clients.clear();
        temp_path.push_back(client);
        temp_queue = client_queue;
        while (!temp_queue.empty() && temp_queue.top() != client) {
            temp_clients.push_back(temp_queue.top());
            temp_queue.pop();
        }
        if (temp_queue.empty()) {
            temp_path.pop_back();
            continue;
        }
        temp_queue.pop();
        for (auto c : temp_clients) {
            temp_queue.push(c);
        }
        res = bruteForceTimesRec(origin, temp_queue, temp_path);
        if (res.first < min_cost) {
            min_cost = res.first;
            min_path = res.second;
        }
        temp_path.pop_back();
    }

    return std::pair<double, std::vector<Vertex<T> *>>(min_cost, min_path);
}

template<class T>
void Graph<T>::nearestNeighbour(const T &origin_info) {
    Vertex<T> *origin = findVertex(origin_info);
    for (auto v : vertexSet)
        v->visited = false;
    origin->visited = true;
    double minWeight;
    Vertex<T> *nextVertex = origin;
    for (unsigned i = 0; i < vertexSet.size() - 1; i++) {
        minWeight = INF;
        auto adj = nextVertex->outgoing;
        for (auto edge : adj) {
            if (!edge->dest->visited && edge->cost < minWeight) {
                minWeight = edge->cost;
                nextVertex = edge->dest;
            }
        }
        nextVertex->visited = true;
        nextVertex->dist = minWeight;
        std::cout << "Dist from " << adj.at(0)->orig->info << " to " << nextVertex->info << ": " << minWeight << "\n";
        nextVertex->path = adj.at(0)->orig;
    }
    origin->path = nextVertex;
    auto adj = nextVertex->outgoing;
    for (auto edge : adj) {
        if (edge->dest == origin) {
            origin->dist = edge->cost;
            break;
        }
    }
}

template<class T>
void Graph<T>::nearestNeighbourTimes(const T &origin_info) {
    Vertex<T> *origin = findVertex(origin_info);
    for (auto v : vertexSet)
        v->visited = false;
    origin->visited = true;
    unsigned min_hour, current_time = start_time;
    double dist;
    Vertex<T> *next_vertex = origin;
    Edge<T> *chosen_edge;
    std::cout << "[Nearest Neighbour Times] Starting search...\n";
    for (unsigned i = 0; i < vertexSet.size() - 1; i++) {
        min_hour = UINT_MAX;
        auto adj = next_vertex->outgoing;
        for (auto edge : adj) {
            if (!edge->dest->visited && edge->dest->hour < min_hour) {
                min_hour = edge->dest->hour;
                next_vertex = edge->dest;
                dist = edge->cost;
                chosen_edge = edge;
            }
        }
        std::cout << "[Nearest Neighbour Times] The time is: " << current_time << ".\n";
        current_time += chosen_edge->cost / velocity;
        std::cout << "[Nearest Neighbour Times] Chose client " << next_vertex->info << ". Arrived at " << current_time
                  << ".\n";
        if (current_time < chosen_edge->dest->hour - early_time) {
            current_time = chosen_edge->dest->hour - early_time;
            std::cout << "[Nearest Neighbour Times] Too early, waiting until " << current_time << ".\n";
        }
        next_vertex->visited_at = current_time;
        current_time += visit_time;
        next_vertex->visited = true;
        next_vertex->dist = dist;
        next_vertex->path = adj.at(0)->orig;
        std::cout << "[Nearest Neighbour Times] Delivered to client " << next_vertex->info << ", leaving at "
                  << current_time << "\n";
    }
    origin->path = next_vertex;
    auto adj = next_vertex->outgoing;
    for (auto edge : adj) {
        if (edge->dest == origin) {
            origin->dist = edge->cost;
            break;
        }
    }
}

//TODO: I think we should change this to a constructor. Or else it just adds on top of the current map...
template<class T>
void Graph<T>::importGraph(std::string vertex_filename, std::string edges_filename, bool lat_lng) {
    latLng = lat_lng;
    std::ifstream vertex_file(vertex_filename);
    std::ifstream edge_file(edges_filename);
    if (vertex_file.fail() || edge_file.fail()) {
        std::cout << "[importGraph] Failed to open files!" << std::endl;
        return;
    }
    unsigned n_vertex = 0;
    char sep;
    unsigned id;
    double lat, lng;
    vertex_file >> n_vertex;
    for (unsigned i = 0; i < n_vertex; i++) {
        vertex_file >> sep;
        vertex_file >> id;
        vertex_file >> sep;
        vertex_file >> lat;
        vertex_file >> sep;
        vertex_file >> lng;
        vertex_file >> sep;
        addVertex(id, lat, lng);
    }

    unsigned n_edges;
    unsigned orig, dest;
    edge_file >> n_edges;
    for (unsigned i = 0; i < n_edges; i++) {
        edge_file >> sep;
        edge_file >> orig;
        edge_file >> sep;
        edge_file >> dest;
        edge_file >> sep;
        addEdge(orig, dest, lat_lng);
    }
}

template<class T>
void Graph<T>::viewGraph() {
    // Instantiate GraphViewer
    GraphViewer gv;
    double factor = 1.0;
    double lat_shift = 0.0;
    double lng_shift = 0.0;
    bool lat_lng = latLng;
    if(lat_lng) factor = LATLNG_FACTOR;

    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));

    unsigned vcounter = 1;
    std::stringstream ss;
    for (auto vertex : getVertexSet()) {
        if(lat_lng){
            lat_shift = vertex->getLat();
            lng_shift = vertex->getLng();
            /*
            if(vertex->getLat() > 0) lat_shift = floor(vertex->getLat());
            else lat_shift = ceil(vertex->getLat());
            if(vertex->getLng() > 0) lng_shift = floor(vertex->getLng());
            else lng_shift = ceil(vertex->getLng());
             */
            lat_lng = false;
        }
        ss << vcounter;
        GraphViewer::Node &node = gv.addNode(vertex->getInfo(), sf::Vector2f((vertex->getLat() - lat_shift) * factor, (vertex->getLng() - lng_shift) * factor));
        node.setLabel(ss.str());
        vcounter++;
    }

    unsigned counter = 0;
    for (auto vertex : getVertexSet()) {
        for (auto edge : vertex->getOutgoing()) {
            gv.addEdge(counter, gv.getNode(vertex->getInfo()), gv.getNode(edge->getDest()->getInfo()),
                       GraphViewer::Edge::DIRECTED);
            counter++;
        }
    }

    // Create window
    gv.createWindow(600, 600);

    // Join viewer thread (blocks till window closed)
    gv.join();
}

template<class T>
void Graph<T>::viewGraphIP(const Graph<T> &igraph) {
    // Instantiate GraphViewer
    GraphViewer gv;
    double factor = 1.0;
    double lat_shift = 0.0;
    double lng_shift = 0.0;
    bool lat_lng = latLng;
    if(lat_lng) factor = LATLNG_FACTOR;

    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));
    bool highlight = false;

    for (auto vertex : getVertexSet()) {
        if(lat_lng){
            lat_shift = vertex->getLat();
            lng_shift = vertex->getLng();
            /*
            if(vertex->getLat() > 0) lat_shift = floor(vertex->getLat());
            else lat_shift = ceil(vertex->getLat());
            if(vertex->getLng() > 0) lng_shift = floor(vertex->getLng());
            else lng_shift = ceil(vertex->getLng());
             */
            lat_lng = false;
        }
        highlight = false;
        GraphViewer::Node &node = gv.addNode(vertex->getInfo(), sf::Vector2f((vertex->getLat() - lat_shift) * factor, (vertex->getLng() - lng_shift) * factor));
        for (auto ivertex : igraph.getVertexSet()) {
            if (ivertex->getInfo() == vertex->getInfo()) {
                highlight = true;
                break;
            }
        }
        if (highlight) {
            node.setColor(GraphViewer::GREEN);
            node.setLabel("INTEREST");
        } else node.setColor(GraphViewer::RED);
    }

    /*
    for(auto ivertex : igraph.getVertexSet()){
        GraphViewer::Node node = gv.addNode(700000 + ivertex->getInfo(), sf::Vector2f(ivertex->getX(), ivertex->getY()));
        node.setColor(GraphViewer::GREEN);
    }
     */

    unsigned counter = 0;
    for (auto vertex : getVertexSet()) {
        for (auto edge : vertex->getOutgoing()) {
            GraphViewer::Edge &road = gv.addEdge(counter, gv.getNode(vertex->getInfo()),
                                                 gv.getNode(edge->getDest()->getInfo()), GraphViewer::Edge::DIRECTED);
            road.setColor(GraphViewer::BLACK);
            counter++;
        }
    }

    for (auto ivertex : igraph.getVertexSet()) {
        for (auto edge : ivertex->getOutgoing()) {
            for (int i = 0; i < (edge->getComplexPath()).size() - 1; i++) {
                for (GraphViewer::Edge *road: gv.getEdges()) {
                    if (road->getFrom()->getId() == (edge->getComplexPath())[i] &&
                        road->getTo()->getId() == (edge->getComplexPath())[i + 1]) {
                        road->setColor(GraphViewer::BLUE);
                        break;
                    }
                }
                //GraphViewer::Edge &road = gv.addEdge(counter, gv.getNode((edge->getComplexPath())[i]), gv.getNode((edge->getComplexPath())[i + 1]), GraphViewer::Edge::DIRECTED);
                //road.setColor(GraphViewer::BLUE);
                //counter++;
            }
        }
    }

    // Create window
    gv.createWindow(600, 600);

    // Join viewer thread (blocks till window closed)
    gv.join();
}

template<class T>
void Graph<T>::viewGraphPath(std::vector<T> path, std::vector<T> interest_points, bool remove_extra_nodes,
                             bool remove_extra_edges, bool show_w, bool show_nid, bool multiple_vans) {
    // Instantiate GraphViewer
    GraphViewer gv;
    std::vector<GraphViewer::Color> color_list = {GraphViewer::BLUE, GraphViewer::GREEN, GraphViewer::RED, GraphViewer::CYAN, GraphViewer::MAGENTA, GraphViewer::LIGHT_GRAY, GraphViewer::YELLOW, GraphViewer::ORANGE, GraphViewer::PINK};
    double factor = 1.0;
    double lat_shift = 0.0;
    double lng_shift = 0.0;
    bool lat_lng = latLng;
    if(lat_lng) factor = LATLNG_FACTOR;
    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));
    bool highlight = false;
    bool special_highlight = false;

    for (auto vertex : getVertexSet()) {
        if(lat_lng){
            lat_shift = vertex->getLat();
            lng_shift = vertex->getLng();
            /*
            if(vertex->getLat() > 0) lat_shift = floor(vertex->getLat());
            else lat_shift = ceil(vertex->getLat());
            if(vertex->getLng() > 0) lng_shift = floor(vertex->getLng());
            else lng_shift = ceil(vertex->getLng());
             */
            lat_lng = false;
        }
        highlight = false;
        special_highlight = false;
        GraphViewer::Node &node = gv.addNode(vertex->getInfo(), sf::Vector2f((vertex->getLat() - lat_shift) * factor, (vertex->getLng() - lng_shift) * factor));
        if (vertex->getInfo() == path[0] || vertex->getInfo() == path[path.size() - 1]) special_highlight = true;
        else {
            for (auto id : interest_points) {
                if (vertex->getInfo() == id) {
                    highlight = true;
                    break;
                }
            }
        }
        if (highlight) {
            node.setColor(GraphViewer::GREEN);
            node.setLabel("INTEREST");
            node.setLabel(std::to_string(node.getId()));
        } else {
            if (special_highlight) {
                node.setLabel("SPECIAL_INTEREST");
                node.setLabel(std::to_string(node.getId()));
            } else {
                node.setColor(GraphViewer::RED);
                if (show_nid) node.setLabel(std::to_string(node.getId()));
            }
        }
    }


    unsigned counter = 0;
    bool special_edge;
    for (auto vertex : getVertexSet()) {
        for (auto edge : vertex->getOutgoing()) {
            special_edge = false;
            GraphViewer::Edge &road = gv.addEdge(counter, gv.getNode(vertex->getInfo()),
                                                 gv.getNode(edge->getDest()->getInfo()), GraphViewer::Edge::DIRECTED);
            for (int i = 0; i < path.size() - 1; i++) {
                if (edge->getOrigin()->getInfo() != 0 && edge->getDest()->getInfo() != 0 && edge->getOrigin()->getInfo() == path[i] && edge->getDest()->getInfo() == path[i + 1]) {
                    special_edge = true;
                    break;
                }
            }
            if (special_edge) {
                if(multiple_vans){
                    GraphViewer::Color color;
                    color = color_list.at(findEdge(road.getFrom()->getId(), road.getTo()->getId())->vanID % color_list.size());
                    road.setColor(color);
                }
                else road.setColor(GraphViewer::BLUE);
                if (show_w) road.setWeight(edge->cost);
            } else road.setColor(GraphViewer::BLACK);
            //road.setWeight(edge->cost);
            counter++;
        }
    }

    /*
    for (int i = 0; i < path.size() - 1; i++) {
        for (GraphViewer::Edge *road: gv.getEdges()) {
            if (road->getFrom()->getId() == path[i] && road->getTo()->getId() == path[i + 1]) {
                road->setColor(GraphViewer::BLUE);
                //road->setLabel(std::to_string(*road->getWeight()));
                break;
            }
        }
    }
     */


    if (remove_extra_nodes) {
        bool ip;
        for (GraphViewer::Node *node : gv.getNodes()) {
            ip = false;
            for (auto id : interest_points) {
                if (id == node->getId()) {
                    ip = true;
                    break;
                }
            }
            //if(!ip) gv.removeNode(node->getId());
            if (!ip) node->setSize(0.0);
        }
    }

    if (remove_extra_edges) {
        bool ie;
        for (GraphViewer::Edge *edge: gv.getEdges()) {
            ie = false;
            for (int i = 0; i < path.size() - 1; i++) {
                if (edge->getFrom()->getId() == path[i] && edge->getTo()->getId() == path[i + 1]) {
                    ie = true;
                    break;
                }
            }
            if (!ie) {
                edge->setThickness(0.0);
            }
        }
    }

    // Create window
    gv.createWindow(600, 600);

    // Join viewer thread (blocks till window closed)
    gv.join();
}


template<class T>
void Graph<T>::viewGraphPathIP(Graph<T> igraph, std::vector<T> path, bool remove_extra_nodes, bool remove_extra_edges,
                               bool show_w, bool show_nid, bool multiple_vans) {
    std::vector<T> interest_points;
    for (auto vertex : igraph.getVertexSet()) {
        interest_points.push_back(vertex->getInfo());
    }

    std::vector<T> real_path;
    for (int i = 0; i < path.size() - 1; i++) {
        Edge<T> *edge = igraph.findEdge(path[i], path[i + 1]);
        for (auto id : edge->complex_path) {
            if (real_path.empty() || real_path[real_path.size() - 1] != id) real_path.push_back(id);
        }
    }

    viewGraphPath(real_path, interest_points, remove_extra_nodes, remove_extra_edges, show_w, show_nid, multiple_vans);
}


template<class T>
bool Graph<T>::checkConectivity(std::vector<T> ids) {
    if (ids.empty()) return false;
    if (ids.size() == 1) return true;
    std::vector<std::vector<T>> sccs = dfsConnectivity();
    std::vector<T> target_scc;
    bool found_target = false;
    for (std::vector<T> scc : sccs) {
        for (auto id : scc) {
            if (id == ids.at(0)) {
                target_scc = scc;
                found_target = true;
                break;
            }
        }
        if (found_target) break;
    }
    if (target_scc.size() < ids.size()) return false;
    bool found_id = false;
    for (int i = 1; i < ids.size(); i++) {
        for (auto id : target_scc) {
            if (id == ids.at(i)) {
                found_id = true;
                break;
            }
        }
        if (!found_id) return false;
    }
    return true;
}

template<class T>
void Graph<T>::printTimes() {
    for (auto vertex : vertexSet) {
        if (vertex->hour > 0) {
            int delta_time = vertex->visited_at - vertex->hour;
            std::cout << "V" << vertex->info << ": tol= " << vertex->tolerance << " | hour= " << vertex->hour
                      << " | visited_at= " << vertex->visited_at
                      << " | delta= " << delta_time << std::endl;
        }
    }
}

template<class T>
std::vector<T> Graph<T>::getOverlapClients(const T &info) {
    std::vector<T> overlap;
    Vertex<T> *vertex = findVertex(info);
    if (vertex == nullptr) return overlap;
    overlap.push_back(info);
    unsigned inf_lim = vertex->hour - early_time;
    unsigned sup_lim = vertex->hour + vertex->max_tolerance;
    for (Vertex<T> *v : vertexSet) {
        if (v->info != vertex->info && v->hour > 0 && v->max_tolerance > 0) {
            //std::cout << "V" << v->info << ": InfTime = " << v->hour - early_time + visit_time << " | SUP_LIM = " << sup_lim << " | SupTime = " << v->hour + v->max_tolerance + visit_time << " | INF_LIM = " << inf_lim << std::endl;
            if ((v->hour >= vertex->hour && v->hour - early_time + visit_time < sup_lim) ||
                (v->hour <= vertex->hour && v->hour + v->max_tolerance + visit_time > inf_lim)) {
                overlap.push_back(v->info);
            }
        }
    }
    return overlap;
}

//Considering an IP Graph (Full Connectivity)
template<class T>
std::vector<T> Graph<T>::getOverlapClientsTravelling(const T &info) {
    std::vector<T> overlap;
    Vertex<T> *vertex = findVertex(info);
    if (vertex == nullptr) return overlap;
    std::vector<Vertex<T>*> res = getOverlapClientsTravelling(vertex);
    for (auto v : res) {
        overlap.push_back(v->info);
    }
    return overlap;
}

template<class T>
std::vector<T> Graph<T>::getPerfectOverlapClientsTravelling(const T &info) {
    std::vector<T> overlap;
    Vertex<T> *vertex = findVertex(info);
    if (vertex == nullptr) return overlap;
    overlap.push_back(info);
    unsigned sup_lim = vertex->hour + vertex->tolerance;
    for (Vertex<T> *v : vertexSet) {
        if (v->info != vertex->info && v->hour > 0 && v->max_tolerance > 0) {
            Edge<T> *edge = findEdge(v, vertex);
            double travelling_time = 0.0;
            if (edge != nullptr) travelling_time = edge->cost / velocity;
            if (v->hour + v->max_tolerance + visit_time + travelling_time < sup_lim) {
                overlap.push_back(v->info);
            }
        }
    }
    return overlap;
}


template<class T>
std::vector<T> Graph<T>::getPerfectOverlapClientsTravelling(std::vector<T> remaining_clients, const T &info) {
    std::vector<T> overlap;
    Vertex<T> *vertex = findVertex(info);
    if (vertex == nullptr) return overlap;
    overlap.push_back(info);
    unsigned sup_lim = vertex->hour + vertex->tolerance;
    for (T vinfo : remaining_clients) {
        Vertex<T> *v = findVertex(vinfo);
        if (v != nullptr && v->info != vertex->info && v->hour > 0 && v->max_tolerance > 0) {
            Edge<T> *edge = findEdge(v, vertex); //TODO: also changed here
            double travelling_time = 0.0;
            if (edge != nullptr) travelling_time = edge->cost / velocity; //TODO: confirm velocity
            if (v->hour + v->max_tolerance + visit_time + travelling_time < sup_lim) {
                overlap.push_back(v->info);
            }
        }
    }
    return overlap;
}

template<class T>
std::vector<Cluster<T>> Graph<T>::getClusters(const T &origin) {
    Vertex<T> *orig = findVertex(origin);
    std::vector<Cluster<T>> result;
    if (orig == nullptr) {
        std::cout << "[Clusters] Invalid origin!\n";
        return result;
    }

    std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareDeparture<T>> cluster_queue;
    //std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareDeparture<T>> temp_queue;
    std::vector<Cluster<T>> temp_clusters;
    std::vector<T> temp;

    for (Vertex<T> *v : vertexSet) {
        temp.clear();
        if (v->info != origin) {
            temp.push_back(v->info);
            Cluster<T> cluster(this, temp);
            cluster_queue.push(cluster);
        }
    }

    Cluster<T> main_cluster = Cluster<T>(this, {});

    while (!cluster_queue.empty()) {
        temp_clusters.clear();
        main_cluster = cluster_queue.top();
        main_cluster.PrintInfo();
        cluster_queue.pop();
        while (!cluster_queue.empty()) {
            Cluster<T> other_cluster = cluster_queue.top();
            if (main_cluster.calculateClusterCost(other_cluster) == 0) {
                main_cluster.mergeCluster(other_cluster);
            } else {
                temp_clusters.push_back(other_cluster);
            }
            cluster_queue.pop();
        }
        for (Cluster<T> c : temp_clusters) cluster_queue.push(c);
        result.push_back(main_cluster);
    }

    return result;
}

/*
template<class T>
std::vector<Cluster<T>> Graph<T>::getClustersCapacity(const T &origin, std::vector<int> capacities) {
    Vertex<T> *orig = findVertex(origin);
    std::vector<Cluster<T>> result;
    if (orig == nullptr) {
        std::cout << "[Clusters] Invalid origin!\n";
        return result;
    }

    std::priority_queue<int> capacities_queue;
    std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareDeparture<T>> cluster_queue;
    //std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareDeparture<T>> temp_queue;
    std::vector<Cluster<T>> temp_clusters;
    std::vector<T> temp;

    for (Vertex<T> *v : vertexSet) {
        temp.clear();
        if (v->info != origin) {
            temp.push_back(v->info);
            Cluster<T> cluster(this, temp);
            cluster_queue.push(cluster);
        }
    }

    Cluster<T> main_cluster = Cluster<T>(this, {});

    while (!cluster_queue.empty()) {
        temp_clusters.clear();
        main_cluster = cluster_queue.top();
        main_cluster.PrintInfo();
        cluster_queue.pop();
        while (!cluster_queue.empty()) {
            Cluster<T> other_cluster = cluster_queue.top();
            if (main_cluster.calculateClusterCost(other_cluster) == 0) {
                main_cluster.mergeCluster(other_cluster);
            } else {
                temp_clusters.push_back(other_cluster);
            }
            cluster_queue.pop();
        }
        for (Cluster<T> c : temp_clusters) cluster_queue.push(c);
        result.push_back(main_cluster);
    }

    return result;
}
*/

template<class T>
class ClusterPairCompare {
public:
    bool operator()(std::pair<std::vector<Cluster<T>>, double> cp1,
                    std::pair<std::vector<Cluster<T>>, double> cp2) {
        return cp1.second > cp2.second;
    }
};

//TODO: Make a Van class and implement it here! (Change everywhere)
template<class T>
std::vector<std::pair<Van, std::vector<Vertex<T> *>>>
Graph<T>::dividingClustersBrute(std::vector<Van> vans, const T &origin) {
    std::cout << "Starting dividingClustersBrute..." << std::endl;
    std::cout << "Getting Clusters..." << std::endl;
    std::vector<Cluster<T>> clusters = getClusters(origin);
    std::vector<Vertex<T> *> vertices;
    std::cout << "Brute Clusters:" << std::endl;
    for(Cluster<T> c : clusters){
        c.PrintInfo();
    }
    unsigned vans_no = vans.size();

    if(vans_no == 0){
        std::vector<std::pair<Van, std::vector<Vertex<T> *>>> empty_result;
        empty_result.push_back(std::pair<Van, std::vector<Vertex<T> *>>(Van(-1, 0, 0), {}));
        return empty_result;
    }
    if(vans_no == 1){
        vertices.clear();
        for(Cluster<T> c : clusters){
            for(T info : c.getInfoSet()){
                vertices.push_back(findVertex(info));
            }
        }
        std::vector<std::pair<Van, std::vector<Vertex<T> *>>> single_result;
        single_result.push_back(std::pair<Van, std::vector<Vertex<T> *>>(vans[0], vertices));
        return single_result;
        //return std::vector<std::pair<Van, std::vector<Vertex<T> *>>>(vans[0], vertices);
    }

    std::vector<Cluster<T>> final_clusters;
    std::vector<std::pair<Van, std::vector<Vertex<T> *>>> result_pairs;

    if (clusters.size() > vans_no) {
        std::cout << "[Brute] More Clusters than Vans..." << std::endl;
        /*
        std::vector<std::vector<double>> cost_matrix(clusters.size() + 1, std::vector<double>(clusters.size() + 1, INF));
        for(int i = 0; i < clusters.size() + 1; i++){
            if(i != 0){
                cost_matrix[i][0] = clusters[i - 1].calculateOriginCost(origin, true);
                cost_matrix[0][i] = clusters[i - 1].calculateOriginCost(origin);
            }
        }

        for(int row = 1; row < clusters.size() + 1; row++){
            for(int col = 1; col < clusters.size() + 1; col++){
                cost_matrix[row][col] = clusters[row - 1].calculateClusterCost(clusters[col - 1]);
            }
        }
         */
        std::vector<double> from_origin_cost(clusters.size(), INF);
        std::vector<double> to_origin_cost(clusters.size(), INF);

        for (int i = 0; i < clusters.size(); i++) {
            from_origin_cost[i] = (clusters[i].calculateOriginCost(origin));
            to_origin_cost[i] = (clusters[i].calculateOriginCost(origin, true));
        }

        int connections = clusters.size() - vans.size();
        std::priority_queue<std::pair<std::vector<Cluster<T>>, double>, std::vector<std::pair<std::vector<Cluster<T>>, double>>, ClusterPairCompare<T>> cluster_connections;
        std::vector<Cluster<T>> connection;


        for (int row = 0; row < clusters.size(); row++) {
            for (int col = 0; col < clusters.size(); col++) {
                if (row != col) {
                    std::vector<double> from_cost = from_origin_cost;
                    std::vector<double> to_cost = to_origin_cost;
                    to_cost[col] = clusters[row].calculateClusterCost(clusters[col]);
                    from_cost[row] = 0.0;
                    double con_cost = std::accumulate(from_cost.begin(), from_cost.end(), 0.0) +
                                      std::accumulate(to_cost.begin(), to_cost.end(), 0.0);
                    connection.clear();
                    connection.push_back(clusters[row]);
                    connection.push_back(clusters[col]);
                    cluster_connections.push(std::pair<std::vector<Cluster<T>>, double>(connection, con_cost));
                }
            }
        }

        std::vector<std::pair<Cluster<T>, Cluster<T> >> cluster_pairs;
        for (int i = 0; i < connections; i++) {
            std::vector<Cluster<T>> cluster_con = cluster_connections.top().first;
            cluster_connections.pop();
            cluster_pairs.push_back(std::pair<Cluster<T>, Cluster<T> >(cluster_con[0], cluster_con[1]));
        }
        for (Cluster<T> cluster : clusters) {
            bool found = false;
            for (std::pair<Cluster<T>, Cluster<T> > cp : cluster_pairs) {
                if (cluster.getInfoSet() == cp.first.getInfoSet() || cluster.getInfoSet() == cp.second.getInfoSet()) {
                    found = true;
                    break;
                }
            }
            if (!found) final_clusters.push_back(cluster);
        }

        /*
        for(auto i = cluster_pairs.begin(); i < cluster_pairs.end() - 1; i++){
            for(auto j = i + 1; j < cluster_pairs.end(); j++){
                std::pair< Cluster<T>, Cluster<T> > cp1 = cluster_pairs[i];
                std::pair< Cluster<T>, Cluster<T> > cp2 = cluster_pairs[j];
                if(cp1.first == cp2.first || cp1.first == cp2.second || cp1.second == cp2.first || cp1.second == cp2.second){
                    Cluster<T> merge_cluster = cp1.first.mergeCluster(cp1.second);
                    merge_cluster.mergeCluster(cp2.first);
                    merge_cluster.mergeCluster(cp2.second);

                }
            }
        }
         */
        std::vector<Cluster<T>> merging_pairs;
        for (int i = 0; i < cluster_pairs.size(); i++) {
            (cluster_pairs[i].first).mergeCluster(cluster_pairs[i].second);
            merging_pairs.push_back(cluster_pairs[i].first);
            //merging_pairs.push_back((cluster_pairs[i].first).mergeCluster(cluster_pairs[i].second));
        }

        while (!merging_pairs.empty()) {
            Cluster<T> c = merging_pairs.back();
            merging_pairs.pop_back();
            bool merged = false;
            for (int i = 0; i < merging_pairs.size(); i++) {
                if (merging_pairs[i].hasCommon(c)) {
                    (merging_pairs[i]).mergeCluster(c);
                    merged = true;
                    break;
                }
            }
            if (!merged) final_clusters.push_back(c);
        }

    } else {
        std::cout << "[Brute] Less Clusters than Vans..." << std::endl;
        final_clusters = clusters;
    }

    std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareSize<T>> final_cluster_queue;
    std::priority_queue<Van, std::vector<Van>, VanCompare> van_queue;
    for (Van van : vans) {
        van_queue.push(van);
    }

    for (Cluster<T> cluster : final_clusters) {
        final_cluster_queue.push(cluster);
    }
    //int vcounter = 0;
    while (!final_cluster_queue.empty()) {
        Cluster<T> cluster = final_cluster_queue.top();
        final_cluster_queue.pop();
        vertices.clear();
        /*
        for (T info : cluster.getInfoSet()) {
            Vertex<T> *v = findVertex(info);
            if (v != nullptr) vertices.push_back(v);
        }
         */
        vertices = cluster.exportVertices();
        //result_pairs.push_back(std::pair<Van, std::vector<Vertex<T> *>>(vans[vcounter], vertices));
        result_pairs.push_back(std::pair<Van, std::vector<Vertex<T> *>>(van_queue.top(), vertices));
        van_queue.pop();
    }

    return result_pairs;
}


template<class T>
std::vector<std::pair<Van, std::vector<Vertex<T> *>>>
Graph<T>::dividingClustersGreedy(std::vector<Van> vans, const T &origin) {
    std::cout << "Starting dividingClustersGreedy..." << std::endl;
    std::cout << "Getting Clusters..." << std::endl;
    std::vector<Cluster<T>> clusters = getClusters(origin);
    std::vector<Vertex<T> *> vertices;
    std::cout << "Greedy Clusters:" << std::endl;
    for(Cluster<T> c : clusters){
        c.PrintInfo();
    }
    unsigned vans_no = vans.size();

    if(vans_no == 0){
        std::vector<std::pair<Van, std::vector<Vertex<T> *>>> empty_result;
        empty_result.push_back(std::pair<Van, std::vector<Vertex<T> *>>(Van(-1, 0, 0), {}));
        return empty_result;
    }
    if(vans_no == 1){
        vertices.clear();
        for(Cluster<T> c : clusters){
            for(T info : c.getInfoSet()){
                vertices.push_back(findVertex(info));
            }
        }
        std::vector<std::pair<Van, std::vector<Vertex<T> *>>> single_result;
        single_result.push_back(std::pair<Van, std::vector<Vertex<T> *>>(vans[0], vertices));
        return single_result;
        //return std::vector<std::pair<Van, std::vector<Vertex<T> *>>>(vans[0], vertices);
    }

    std::vector<Cluster<T>> final_clusters;
    std::vector<std::pair<Van, std::vector<Vertex<T> *>>> result_pairs;

    if (clusters.size() > vans_no) {
        std::cout << "[Greedy] More Clusters than Vans..." << std::endl;
        std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareArrival<T>> arrival_cluster_queue;
        std::vector<Cluster<T>> van_clusters(vans_no, Cluster<T>(this, {}));
        for(int i = 0; i < clusters.size(); i++){
            arrival_cluster_queue.push(clusters.at(i));
        }
        int vcounter = 0;
        while(!arrival_cluster_queue.empty()){
            //van_clusters[vcounter] = (van_clusters[vcounter]).mergeCluster(arrival_cluster_queue.top());
            (van_clusters[vcounter]).mergeCluster(arrival_cluster_queue.top());
            van_clusters[vcounter].PrintInfo();
            arrival_cluster_queue.pop();
            vcounter++;
            if(vcounter >= vans_no) vcounter = 0;
        }
        final_clusters = van_clusters;
    } else {
        std::cout << "[Greedy] Less Clusters than Vans..." << std::endl;
        final_clusters = clusters;
    }

    std::priority_queue<Cluster<T>, std::vector<Cluster<T>>, ClusterCompareSize<T>> final_cluster_queue;
    std::priority_queue<Van, std::vector<Van>, VanCompare> van_queue;
    for (Van van : vans) {
        van_queue.push(van);
    }

    for (Cluster<T> cluster : final_clusters) {
        final_cluster_queue.push(cluster);
    }
    //int vcounter = 0;
    while (!final_cluster_queue.empty()) {
        Cluster<T> cluster = final_cluster_queue.top();
        final_cluster_queue.pop();
        vertices.clear();
        /*
        for (T info : cluster.getInfoSet()) {
            Vertex<T> *v = findVertex(info);
            if (v != nullptr) vertices.push_back(v);
        }
         */
        vertices = cluster.exportVertices();
        //result_pairs.push_back(std::pair<Van, std::vector<Vertex<T> *>>(vans[vcounter], vertices));
        result_pairs.push_back(std::pair<Van, std::vector<Vertex<T> *>>(van_queue.top(), vertices));
        van_queue.pop();
    }

    return result_pairs;
}


template<class T>
std::vector<Vertex<T> *> Graph<T>::getOverlapClientsTravelling(Vertex<T> *info) {
    std::vector<Vertex<T> *> overlap;
    if (info == nullptr) return overlap;
    overlap.push_back(info);
    unsigned inf_lim = info->hour - early_time;
    unsigned sup_lim = info->hour + info->max_tolerance;
    for (Vertex<T> *v : vertexSet) {
        if (v->info != info->info && v->hour > 0 && v->max_tolerance > 0) {
            Edge<T> *edge = findEdge(v, info);
            double travelling_time = 0.0;
            if (edge != nullptr) travelling_time = edge->cost / velocity;
            if ((v->hour >= info->hour && v->hour - early_time + visit_time + travelling_time < sup_lim) ||
                (v->hour <= info->hour && v->hour + v->max_tolerance + visit_time + travelling_time > inf_lim)) {
                overlap.push_back(v);
            }
        }
    }
    return overlap;
}


template<class T>
double Graph<T>::costFunctionStep(T og, std::vector<T> path, T new_element) {
    if (weight < 0 || weight > 1) return -1;

    unsigned current_time = start_time;
    double delay = 0;
    double total_delay = 0;

    double average; //f
    double deviation; //g

    T first;
    T second;
    Edge<T> *edge;
    Vertex<T> *current_v;

    Vertex<T> *origin = findVertex(og);
    Vertex<T> *new_vertex = findVertex(new_element);
    if (origin == nullptr || new_vertex == nullptr) return -1;

    for (int i = 0; i < path.size(); i++) {
        if (i == 0) first = origin->info;
        else first = path[i - 1];
        second = path[i];
        edge = findEdge(first, second);
        current_v = findVertex(second);
        if (edge == nullptr || current_v == nullptr) return -1;
        delay = (current_time + edge->cost / velocity) - current_v->hour;
        if (delay < current_v->tolerance) delay = 0;
        if (delay > current_v->max_tolerance) return INF;
        total_delay += delay;
        if (current_time + edge->cost / velocity < current_v->hour - early_time)
            current_time = current_v->hour - early_time + visit_time;
        else current_time += edge->cost / velocity + visit_time;
    }
    if (path.size() > 0) edge = findEdge(path[path.size() - 1], new_vertex->info);
    else edge = findEdge(origin, new_vertex);
    if (edge == nullptr) return -1;
    delay = (current_time + edge->cost / velocity) - new_vertex->hour;
    if (delay < new_vertex->tolerance) delay = 0;
    if (delay > new_vertex->max_tolerance) return INF;
    total_delay += delay;
    double vertex_no = path.size() + 1;
    average = total_delay / vertex_no;
    deviation = sqrt((1 / vertex_no) * pow(delay - average, 2));

    double cost = weight * average + (1 - weight) * deviation;
    return cost;
}


template<class T>
double Graph<T>::costFunctionTotal(T og, std::vector<T> path) {

    Vertex<T> *origin = findVertex(og);
    std::vector<Vertex<T> *> ptr_path;
    for (auto v : path) {
        ptr_path.push_back(findVertex(v));
    }
    return costFunctionTotal(origin, ptr_path);
}

template<class T>
double Graph<T>::costFunctionTotal(Vertex<T> *origin, std::vector<Vertex<T> *> path) {
    if (weight < 0 || weight > 1) return -1;
    if (path.size() == 0) return -1;

    std::vector<unsigned> delays;

    unsigned current_time = start_time;
    double delay = 0;
    double total_delay = 0;

    double average; //f
    double deviation; //g

    Vertex<T> *first;
    Vertex<T> *second;
    Edge<T> *edge;
    Vertex<T> *current_v;

    if (origin == nullptr) return -1;

    for (int i = 0; i < path.size(); i++) {
        if (i == 0) first = origin;
        else first = path[i - 1];
        second = path[i];
        edge = findEdge(first, second);
        current_v = second;
        if (edge == nullptr || current_v == nullptr) return -1;
        delay = (current_time + edge->cost / velocity) - current_v->hour;
        if (delay < current_v->tolerance) delay = 0;
        if (delay > current_v->max_tolerance) return INF;
        delays.push_back(delay);
        total_delay += delay;
        if (current_time + edge->cost / velocity < current_v->hour - early_time)
            current_time = current_v->hour - early_time + visit_time;
        else current_time += edge->cost / velocity + visit_time;
    }

    double vertex_no = path.size();
    average = total_delay / vertex_no;

    double deviation_sum = 0;
    for (int i = 0; i < delays.size(); i++) {
        deviation_sum += pow(delays[i] - average, 2);
    }
    deviation = sqrt((1 / vertex_no) * deviation_sum);

    double cost = weight * average + (1 - weight) * deviation;
    return cost;
}

template<class T>
void Graph<T>::followVansPath(std::vector<std::pair<Van, std::vector<Vertex<T> *>>> paths_dist, const T &og){
    Vertex<T> * origin = findVertex(og);
    if (origin == nullptr) {
        std::cout << "[Van Paths] Invalid origin!\n";
        return;
    }
    for(int i = 0; i < paths_dist.size(); i++){
        double current_time = start_time;
        double visit_time = paths_dist[i].first.getVisitTime();
        double travelling_time = 0;
        int van_id = paths_dist[i].first.getID();
        Vertex<T> * current_client;
        Edge<T> * road;

        //"[Van Paths] The time is: "
        /*
        std::cout << "[Van Paths] Van " << van_id << "-> The time is: " << current_time << ".\n";
        std::cout << "[Van Paths] Van " << van_id << "-> Chose client " << next_vertex->info << ". Arrived at " << current_time << ".\n";
        std::cout << "[Van Paths] Van " << van_id << "-> Too early, waiting until " << current_time << ".\n";
        std::cout << "[Van Paths] Van " << van_id << "-> Delivered to client " << next_vertex->info << ", leaving at " << current_time << "\n";
        */
        for(int j = 0; j < paths_dist[i].second.size(); j++){
            std::cout << "[Van Paths] Van " << van_id << "-> The time is: " << current_time << ".\n";
            current_client = paths_dist[i].second.at(j);
            if(j == 0) road = findEdge(origin, current_client);
            else road = findEdge(paths_dist[i].second.at(j-1), current_client);
            if(road == nullptr){
                std::cout << "[Van Paths] Invalid Connection!\n";
                return;
            }
            travelling_time = road->cost / velocity;
            current_time += travelling_time;
            std::cout << "[Van Paths] Van " << van_id << "-> Chose client " << current_client->info << ". Arrived at " << current_time << ".\n";
            if(current_time < current_client->hour - early_time){
                current_time = current_client->hour - early_time;
                std::cout << "[Van Paths] Van " << van_id << "-> Too early, waiting until " << current_time << ".\n";
            }
            double delay = current_time - current_client->hour + current_client->tolerance;
            if(delay < 0) delay = 0;
            if(delay == 0) std::cout << "[Van Paths] Van " << van_id << "-> Arrived perfectly on time" << ".\n";
            else std::cout << "[Van Paths] Van " << van_id << "-> Arrived " << delay << " minutes late" << ".\n";
            current_client->visited = true;
            current_client->visited_at = current_time;
            road->vanID = van_id;
            current_time += visit_time;
            std::cout << "[Van Paths] Van " << van_id << "-> Delivered to client " << current_client->info << ", leaving at " << current_time << "\n";
        }
        road = findEdge(current_client, origin);
        if(road == nullptr){
            std::cout << "[Van Paths] Invalid Connection!\n";
            return;
        }
        road->vanID = van_id;
        travelling_time = road->cost / velocity;
        current_time += travelling_time;
        std::cout << "[Van Paths] Van " << van_id << "-> Returning to the Origin. Arrived, finally, at " << current_time << "\n";

    }
}

template<class T>
std::vector<T> Graph<T>::vanPathtoViewPath(std::vector<std::pair<Van, std::vector<Vertex<T > *>>> van_path, const T &og){
    std::vector<T> result;
    for(int i = 0; i < van_path.size(); i++){
        result.push_back(og);
        for(int j = 0; j < van_path[i].second.size(); j++){
            result.push_back(van_path[i].second[j]->info);
        }
        result.push_back(og);
    }
    return result;
}

#endif /* PROJ_GRAPH_H_ */
