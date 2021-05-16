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

    bool visited;  // for path finding
    //Edge<T> *path; // for path finding
    Vertex<T> *path = NULL;
    double dist;   // for path finding
    int queueIndex = 0; // required by MutablePriorityQueue

    Vertex(T in);

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

    friend class Graph<T>;

    friend class MutablePriorityQueue<Vertex<T>>;
};

template<class T>
Vertex<T>::Vertex(T in): info(in) {}

template<class T>
Vertex<T>::Vertex(T in, double x, double y): info(in), x(x), y(y) {}

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
    double capacity;
    double cost;
    double flow;

    Edge(Vertex<T> *o, Vertex<T> *d, double capacity, double cost = 0, double flow = 0);

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
Edge<T>::Edge(Vertex<T> *o, Vertex<T> *d, double capacity, double cost, double flow):
        orig(o), dest(d), capacity(capacity), cost(cost), flow(flow) {}

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

public:
    std::vector<Vertex<T> *> nearestNeighbour(const T &origin_info);

    void dijkstraShortestPath(const T &origin);

    void heldKarp(const T &origin);

    std::vector<std::vector<T>> dfsConnectivity() const;

    bool checkConectivity(std::vector<T> ids);

    std::vector<T> dfs() const;

    void dfsVisit(Vertex<T> *v, std::vector<T> &res, bool reverse = false) const;

    Graph<T> generateInterestPointsGraph(std::vector<T> important_points);

    Vertex<T> *findVertex(const T &inf) const;

    std::vector<Vertex<T> *> getVertexSet() const;

    Vertex<T> *addVertex(const T &in);

    Vertex<T> *addVertex(const T &in, const double &x, const double &y);

    Edge<T> *addEdge(const T &sourc, const T &dest, double capacity, double cost, double flow = 0);

    Edge<T> *addEdge(const T &sourc, const T &dest, double cost, std::vector<T> complex_path);

    Edge<T> *addEdge(const T &sourc, const T &dest, double cost);

    Edge<T> *addEdge(const T &sourc, const T &dest);

    std::vector<T> getPath(const T &origin, const T &dest) const;

    void viewGraph();

    void viewGraphPath(Graph<T> igraph);

    void importGraph(std::string vertex_filename, std::string edges_filename);
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
std::vector<Vertex<T> *> Graph<T>::getVertexSet() const {
    return vertexSet;
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
std::vector<T> Graph<T>::getPath(const T &origin, const T &dest) const {
    std::vector<T> res;

    Vertex<T> *v = this->findVertex(dest);
    if (v == nullptr || v->dist == INF) return res;

    while (true) {
        res.push_back(v->info);
        if (v->info == origin) break;
        v = v->path;
    }

    std::reverse(res.begin(), res.end());

    return res;
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
        std::cout << std::endl << "I(" << i << ") -> currentInfo = " << current_info << std::endl;
        for (int j = 0; j < important_points.size(); j++) {
            if (important_points[j] != current_info) {
                std::vector<T> path = this->getPath(current_info, important_points[j]);
                double dist = this->findVertex(important_points[j])->getDist();
                std::cout << "J(" << j << ") -> Info = " << important_points[j] << std::endl;
                std::cout << "J(" << j << ") -> dist = " << dist << std::endl;
                std::cout << "J(" << j << ") -> path = ";
                for (int a = 0; a < path.size(); a++) std::cout << path[a] << ", ";
                std::cout << std::endl;
                result.addEdge(current_info, important_points[j], dist, path);
            }
        }
    }
    return result;
}

template<class T>
void Graph<T>::heldKarp(const T &origin) {

}

template<class T>
std::vector<Vertex<T> *> Graph<T>::nearestNeighbour(const T &origin_info) { //TODO: Test!
    Vertex<T> *origin = findVertex(origin_info);
    for (auto v : vertexSet)
        v->visited = false;
    std::vector<Vertex<T> *> path;
    origin->visited = true;
    path.push_back(origin);
    double minWeight;
    Vertex<T> *nextVertex = origin;
    do {
        minWeight = INF;
        auto adj = nextVertex->outgoing;
        for (auto edge : adj) {
            if (!edge->dest->visited && edge->weight < minWeight) {
                minWeight = edge->cost;
                nextVertex = edge->dest;
            }
        }
        nextVertex->visited = true;
        path.push_back(nextVertex);
    } while (minWeight != UINT_MAX);
    path.push_back(origin);
    return path;
}

template<class T>
void Graph<T>::importGraph(std::string vertex_filename, std::string edges_filename) {
    std::ifstream vertex_file(vertex_filename);
    std::ifstream edge_file(edges_filename);
    if (vertex_file.fail() || edge_file.fail()) {
        std::cout << "Failed to open files!" << std::endl;
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
void Graph<T>::viewGraphPath(Graph<T> igraph) {
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
            std::cout << "HIGHLIGHT" << std::endl;
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

#endif /* PROJ_GRAPH_H_ */
