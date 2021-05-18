#ifndef PROJ_GRAPH_H_
#define PROJ_GRAPH_H_

#include <vector>
#include <queue>
#include <limits>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <graphviewer.h>
#include <sstream>
#include <string>
#include "MutablePriorityQueue.h"
#include "Utils.h"


constexpr auto INF = std::numeric_limits<double>::max();

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
    T info;
    double x;
    double y;
    std::vector<Edge<T> *> outgoing;
    std::vector<Edge<T> *> incoming;
    std::vector<Edge<T> *> adj;

    bool visited = false;  // for path finding
    //Edge<T> *path; // for path finding
    Vertex<T> *path = NULL;
    double dist = INF;   // for path finding
    int queueIndex = 0; // required by MutablePriorityQueue

    unsigned hour;
    unsigned tolerance;
    unsigned max_tolerance;
    unsigned quantity;
    unsigned visited_at;

    explicit Vertex(T in);

    Vertex(T in, double x, double y, unsigned hour, unsigned tolerance, unsigned max_tolerance);

    Vertex(T in, double x, double y);

    void addEdge(Edge<T> *e);

    bool operator<(Vertex<T> &vertex) const; // required by MutablePriorityQueue

public:
    T getInfo() const;

    double getDist() { return dist; }

    std::vector<Edge<T> *> getIncoming() const;

    std::vector<Edge<T> *> getOutgoing() const;

    double getX() const { return x; }

    double getY() const { return y; }

    void setX(double x) { this->x = x; }

    void setY(double y) { this->y = y; }

    void setTimes(unsigned hour, unsigned tolerance, unsigned max_tolerance) {this->hour = hour; this->tolerance = tolerance; this->max_tolerance = max_tolerance;}

    friend class Graph<T>;

    friend class MutablePriorityQueue<Vertex<T>>;
};

template<class T>
Vertex<T>::Vertex(T in): info(in), x(0), y(0) {}

template<class T>
Vertex<T>::Vertex(T in, double x, double y): info(in), x(x), y(y) {}

template<class T>
Vertex<T>::Vertex(T in, double x, double y, unsigned hour, unsigned tolerance, unsigned max_tolerance): info(in), x(x), y(y), hour(hour), tolerance(tolerance), max_tolerance(max_tolerance) {}

template<class T>
void Vertex<T>::addEdge(Edge<T> *e) {
    adj.push_back(e);
    outgoing.push_back(e);
    e->dest->incoming.push_back(e);
    e->dest->adj.push_back(e);
}

template<class T>
bool Vertex<T>::operator<(Vertex<T> &vertex) const {
    return this->dist < vertex.dist;
}

template<class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template<class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

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
    Vertex<T> *orig;
    Vertex<T> *dest;
    std::vector<T> complex_path;
    double cost;

public:
    Edge(Vertex<T> *o, Vertex<T> *d, double cost, std::vector<T> complex_path);

    Edge(Vertex<T> *o, Vertex<T> *d, double cost);

    Edge(Vertex<T> *o, Vertex<T> *d);

    friend class Graph<T>;

    friend class Vertex<T>;

    double getFlow() const;

    //std::vector<Edge<T> *> getComplexPath();
    std::vector<T> getComplexPath();

    void setComplexPath(std::vector<Edge<T> *> complex_path);

    Vertex<T> *getOrigin() { return orig; }

    Vertex<T> *getDest() { return dest; }

    double getCost() const { return cost; }

    void setCost(double cost) { this->cost = cost; }
};

template<class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double cost, std::vector<T> complex_path): orig(o), dest(d), cost(cost) {
    this->complex_path = complex_path;
}

template<class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double cost):
        orig(o), dest(d), cost(cost) {}

template<class T>
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d):
        orig(o), dest(d) {
    cost = calculateDist(this->orig->getX(), this->orig->getY(), this->dest->getX(), this->dest->getY());
}

template<class T>
double Edge<T>::getFlow() const {
    return this->flow;
}

/*
template<class T>
std::vector<Edge<T> *> Edge<T>::getComplexPath() {
    return this->complex_path;
}
 */

template<class T>
std::vector<T> Edge<T>::getComplexPath() {
    return this->complex_path;
}

template<class T>
void Edge<T>::setComplexPath(std::vector<Edge<T> *> complex_path) {
    this->complex_path = complex_path;
}

/*
 * ================================================================================================
 * Class Graph
 * ================================================================================================
 */

template<class T>
class Graph {
    std::vector<Vertex<T> *> vertexSet;
    bool ipGraph = false; //Interest points graph
    std::vector<std::vector<double>> distance; //Initialize only for IP Graphs
    unsigned start_time = 0;
public:
    unsigned int getStartTime() const;

    void setStartTime(unsigned int startTime);

    unsigned int getVisitTime() const;

    void setVisitTime(unsigned int visitTime);

    unsigned int getEarlyTime() const;

    void setEarlyTime(unsigned int earlyTime);

private:
    unsigned visit_time = 5;
    unsigned early_time = 10;
    unsigned velocity = 50;

public:
    void nearestNeighbour(const T &origin_info);

    void nearestNeighbourTimes(const T &origin_info);

    void dijkstraShortestPath(const T &origin);

    void heldKarp(const T &origin);

    std::vector<std::vector<T>> dfsConnectivity() const;

    bool checkConectivity(std::vector<T> ids);

    std::vector<T> dfs() const;

    void dfsVisit(Vertex<T> *v, std::vector<T> &res, bool reverse = false) const;

    Graph<T> generateInterestPointsGraph(std::vector<T> important_points);

    Vertex<T> *findVertex(const T &inf) const;

    Edge<T> *findEdge(const T &og, const T &dest) const;

    std::vector<Vertex<T> *> getVertexSet() const;

    Vertex<T> *addVertex(const T &in);

    Vertex<T> *addVertex(const T &in, const double &x, const double &y);

    Vertex<T> *addVertex(const T &in, const double &x, const double &y, const unsigned &hour, const unsigned &tolerance, const unsigned &max_tolerance);

    Edge<T> *addEdge(const T &sourc, const T &dest, double capacity, double cost, double flow = 0);

    Edge<T> *addEdge(const T &sourc, const T &dest, double cost, std::vector<T> complex_path);

    Edge<T> *addEdge(const T &sourc, const T &dest, double cost);

    Edge<T> *addEdge(const T &sourc, const T &dest);

    std::pair<std::vector<T> , double> getPath(const T &origin, const T &dest) const;

    unsigned getVelocity() {return this->velocity;}

    void setVelocity(unsigned velocity) {this->velocity = velocity;}

    void viewGraph();

    void viewGraphIP(Graph<T> igraph);

    void viewGraphPath(std::vector<T> path, std::vector<T> interest_points, bool remove_extra_nodes = false,
                       bool remove_extra_edges = false, bool show_w=false, bool show_nid=false);

    void viewGraphPathIP(Graph<T> igraph, std::vector<T> path, bool remove_extra_nodes = false,
                         bool remove_extra_edges = false, bool show_w=false, bool show_nid=false);

    void importGraph(std::string vertex_filename, std::string edges_filename);

    void printTimes();
};

template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in);
    vertexSet.push_back(v);
    return v;
}

template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in, const double &x, const double &y) {
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in, x, y);
    vertexSet.push_back(v);
    return v;
}

template<class T>
Vertex<T> *Graph<T>::addVertex(const T &in, const double &x, const double &y, const unsigned &hour, const unsigned &tolerance, const unsigned &max_tolerance){
    Vertex<T> *v = findVertex(in);
    if (v != nullptr)
        return v;
    v = new Vertex<T>(in, x, y, hour, tolerance, max_tolerance);
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
Edge<T> *Graph<T>::addEdge(const T &sourc, const T &dest) {
    auto s = findVertex(sourc);
    auto d = findVertex(dest);
    if (s == nullptr || d == nullptr)
        return nullptr;
    Edge<T> *e = new Edge<T>(s, d);
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
 * *************** DFS  ************
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
 * *************** Shortest Path Problem  ***********
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
std::pair<std::vector<T> , double> Graph<T>::getPath(const T &origin, const T &dest) const {
    std::vector<T> path;
    double total_dist = 0;

    Vertex<T> *v = this->findVertex(dest);
    if (v == nullptr || v->dist == INF) return std::pair<std::vector<T> , double>(path, total_dist);

    while (true) {
        path.push_back(v->info);
        total_dist += v->dist;
        v = v->path;
        if (v == nullptr) {
            return std::pair<std::vector<T> , double>(path, total_dist);
        }
        if (v->info == origin) {
            path.push_back(v->info);
            break;
        }
    }

    std::reverse(path.begin(), path.end());

    return std::pair<std::vector<T> , double>(path, total_dist);
}

template<class T>
Graph<T> Graph<T>::generateInterestPointsGraph(std::vector<T> important_points) {
    Graph<T> result;
    std::cout << "Checking Con..." << std::endl;
    if (!checkConectivity(important_points)) {
        std::cout << "Conectivity Invalid!!!" << std::endl;
        return result;
    }
    std::cout << "Con Good!" << std::endl;

    for (int i = 0; i < important_points.size(); i++) {
        Vertex<T> *v = findVertex(important_points[i]);
        result.addVertex(important_points[i], v->getX(), v->getY());
    }
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

template<class T>
void Graph<T>::heldKarp(const T &origin) {
    //Initialize vertices vector with origin in last place
    Vertex<T> *orig = findVertex(origin);
    if (orig == nullptr) {
        std::cout << "[Held-Karp] Invalid origin!\n";
        return;
    }
    std::vector<Vertex<T> *> vertices = vertexSet;
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
    for (unsigned visited = 1; visited < (1 << vertices.size() - 1); visited++) {
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
    unsigned final_dist, second_to_last;
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
        std::cout << "[Nearest Neighbour Times] Choosed client " << next_vertex->info << ". Arrived at" << current_time << ".\n";
        if (current_time < chosen_edge->dest->hour - early_time) {
            current_time = chosen_edge->dest->hour - early_time;
            std::cout << "[Nearest Neighbour Times] Too early, waiting until " << current_time << ".\n";
        }
        current_time += visit_time;
        next_vertex->visited = true;
        next_vertex->dist = dist;
        next_vertex->visited_at = current_time;
        next_vertex->path = adj.at(0)->orig;
        std::cout << "[Nearest Neighbour Times] Delivered to client " << next_vertex->info << ", leaving at " << current_time << "\n";
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

template<class T>
void Graph<T>::importGraph(std::string vertex_filename, std::string edges_filename) {
    std::ifstream vertex_file(vertex_filename);
    std::ifstream edge_file(edges_filename);
    if (vertex_file.fail() || edge_file.fail()) {
        std::cout << "[importGraph] Failed to open files!" << std::endl;
        return;
    }
    unsigned n_vertex = 0;
    char sep;
    unsigned id;
    double x, y;
    vertex_file >> n_vertex;
    for (unsigned i = 0; i < n_vertex; i++) {
        vertex_file >> sep;
        vertex_file >> id;
        vertex_file >> sep;
        vertex_file >> x;
        vertex_file >> sep;
        vertex_file >> y;
        vertex_file >> sep;
        addVertex(id, x, y);
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
        addEdge(orig, dest);
    }
}

template<class T>
void Graph<T>::viewGraph() {
    // Instantiate GraphViewer
    GraphViewer gv;

    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));

    unsigned vcounter = 1;
    std::stringstream ss;
    for (auto vertex : getVertexSet()) {
        ss << vcounter;
        GraphViewer::Node &node = gv.addNode(vertex->getInfo(), sf::Vector2f(vertex->getX(), vertex->getY()));
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
void Graph<T>::viewGraphIP(Graph<T> igraph) {
    // Instantiate GraphViewer
    GraphViewer gv;

    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));
    bool highlight = false;

    for (auto vertex : getVertexSet()) {
        highlight = false;
        GraphViewer::Node &node = gv.addNode(vertex->getInfo(), sf::Vector2f(vertex->getX(), vertex->getY()));
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
                             bool remove_extra_edges, bool show_w, bool show_nid) {
    // Instantiate GraphViewer
    GraphViewer gv;
    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));
    bool highlight = false;
    bool special_highlight = false;

    for (auto vertex : getVertexSet()) {
        highlight = false;
        special_highlight = false;
        GraphViewer::Node &node = gv.addNode(vertex->getInfo(), sf::Vector2f(vertex->getX(), vertex->getY()));
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
            } else{
                node.setColor(GraphViewer::RED);
                if(show_nid) node.setLabel(std::to_string(node.getId()));
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
                if (edge->getOrigin()->getInfo() == path[i] && edge->getDest()->getInfo() == path[i + 1]){
                    special_edge = true;
                    break;
                }
            }
            if(special_edge){
                road.setColor(GraphViewer::BLUE);
                if(show_w) road.setWeight(edge->cost);
            }
            else road.setColor(GraphViewer::BLACK);
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
            for(int i = 0; i < path.size() - 1; i++){
                if (edge->getFrom()->getId() == path[i] && edge->getTo()->getId() == path[i + 1]) {
                    ie = true;
                    break;
                }
            }
            if(!ie) {
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
void Graph<T>::viewGraphPathIP(Graph<T> igraph, std::vector<T> path, bool remove_extra_nodes, bool remove_extra_edges, bool show_w, bool show_nid) {
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

    viewGraphPath(real_path, interest_points, remove_extra_nodes, remove_extra_edges, show_w, show_nid);
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
void Graph<T>::printTimes(){
    for(auto vertex : vertexSet){
        if(vertex->hour > 0){
            unsigned delta_time = vertex->visited_at - vertex->hour;
            std::cout << "V" << vertex->info << ": tol= " << vertex->tolerance << " | max_tol= " << vertex->max_tolerance << " | hour= " << vertex->hour << " | visited_at= " << vertex->visited_at << " | delta= " << delta_time << std::endl;
        }
    }
}


#endif /* PROJ_GRAPH_H_ */
