#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <queue>
#include <limits>
#include <iostream>
#include <algorithm>
#include "MutablePriorityQueue.h"

constexpr auto INF = std::numeric_limits<double>::max();

template <class T> class Vertex;
template <class T> class Edge;
template <class T> class Graph;

/*
 * ================================================================================================
 * Class Vertex
 * ================================================================================================
 */

template <class T>
class Vertex {
	T info;
	std::vector<Edge<T> *> outgoing;
    std::vector<Edge<T> *> incoming;
    std::vector<Edge<T> *> adj;

	bool visited;  // for path finding
	//Edge<T> *path; // for path finding
    Vertex<T> *path = NULL;
	double dist;   // for path finding
	int queueIndex = 0; // required by MutablePriorityQueue

	Vertex(T in);
	void addEdge(Edge<T> *e);
	bool operator<(Vertex<T> & vertex) const; // required by MutablePriorityQueue

public:
	T getInfo() const;
	double getDist() {return dist;}
    std::vector<Edge<T> *> getIncoming() const;
    std::vector<Edge<T> *> getOutgoing() const;
	friend class Graph<T>;
	friend class MutablePriorityQueue<Vertex<T>>;
};

/*
 * ================================================================================================
 * Class Edge
 * ================================================================================================
 */

template <class T>
class Edge {
	Vertex<T> * orig;
	Vertex<T> * dest;
    std::vector<T> complex_path;
	double capacity;
	double cost;
	double flow;

	Edge(Vertex<T> *o, Vertex<T> *d, double capacity, double cost=0, double flow=0);

public:
    Edge(Vertex<T> *o, Vertex<T> *d, double cost, std::vector<T> complex_path);
	friend class Graph<T>;
	friend class Vertex<T>;
	double getFlow() const;
    std::vector<Edge<T> *> getComplexPath();
    void setComplexPath(std::vector<Edge<T> *> complex_path);
    Vertex<T> * getOrigin() {return orig;}
    Vertex<T> * getDest() {return dest;}
    double getCost() const {return cost;}
};

/*
 * ================================================================================================
 * Class Graph
 * ================================================================================================
 */

template <class T>
class Graph {
    std::vector<Vertex<T> *> vertexSet;

    void dijkstraResidualShortestPath(Vertex<T> *s);
	void bellmanFordShortestPath(Vertex<T> *s);
	bool relax(Vertex<T> *v, Vertex<T> *w, Edge<T> *e, double residual, double cost);

	void resetFlows();
	bool findAugmentationPath(Vertex<T> *s, Vertex<T> *t);
	void testAndVisit(std::queue< Vertex<T>*> &q, Edge<T> *e, Vertex<T> *w, double residual);
	double findMinResidualAlongPath(Vertex<T> *s, Vertex<T> *t);
	void augmentFlowAlongPath(Vertex<T> *s, Vertex<T> *t, double flow);

public:
    std::vector<Vertex<T> *> nearestNeighbour(const T &origin_info);
    void dijkstraShortestPath(const T &origin);
    std::vector<T> dfs() const;
    void dfsVisit(Vertex<T> *v, std::vector<T> & res) const;
    Graph <T> generateInterestPointsGraph(std::vector<T> important_points);
	Vertex<T>* findVertex(const T &inf) const;
    std::vector<Vertex<T> *> getVertexSet() const;
	Vertex<T> *addVertex(const T &in);
	Edge<T> *addEdge(const T &sourc, const T &dest, double capacity, double cost, double flow=0);
    Edge<T> *addEdge(const T &sourc, const T &dest, double cost, std::vector<T> complex_path);
	double getFlow(const T &sourc, const T &dest) const ;
	void fordFulkerson(T source, T target);
	double minCostFlow(T source, T target, double flow);
    std::vector<T> getPath(const T &origin, const T &dest) const;
};

#endif /* GRAPH_H_ */
