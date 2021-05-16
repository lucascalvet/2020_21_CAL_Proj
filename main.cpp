#include <iostream>
#include <string>
#include "Graph.h"
#include <math.h>
using namespace std;

void addPair(vector<int> pairs, int n1, int n2){
    int item = n1 * 100 + n2;
    pairs.push_back(item);
}

bool isPair(vector<int> pairs, int n1, int n2){
    int item = n1 * 100 + n2;
    for(int i = 0; i < pairs.size(); i++){
        if(pairs[i] == item) return true;
    }
    return false;
}

template <class T>
void VertexPrintInfo(Vertex<T> * v) {cout << "V" << v->getInfo();}
template <class T>
void EdgePrintInfo(Edge<T> * e) {VertexPrintInfo(e->getOrigin()); cout << "--" << e->getCost() << "-->"; VertexPrintInfo(e->getDest()); cout << endl;}
template <class T>
void GraphPrintInfo(Graph<T> g) {
    for(Vertex<T>* vertex : g.getVertexSet()){
        for(Edge<T>* edge : vertex->getOutgoing()){
            EdgePrintInfo(edge);
        }
        /*
        for(Edge<T>* edge : vertex->getIncoming()){
            EdgePrintInfo(edge);
        }
         */
    }
}

/*
template <class T>
Graph<T> dijkstraInterestPoints(Graph<T> complete_graph, vector<T> important_points){
    Graph<T> result;
    for(int i = 0; i < important_points.size(); i++){
        result.addVertex(important_points[i]);
    }
    for(int i = 0; i < important_points.size(); i++){
        T current_info = important_points[i];
        //Graph<T> copy = complete_graph;
        complete_graph.dijkstraShortestPath(current_info);
        //copy.dijkstraShortestPath(current_info);
        cout << endl << "I(" << i << ") -> currentInfo = " << current_info << endl;
        for(int j = 0; j < important_points.size(); j++){
            if(important_points[j] != current_info){
                vector<T> path = complete_graph.getPath(current_info, important_points[j]);
                double dist = complete_graph.findVertex(j)->getDist();
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
 */

int main(){
    cout << "Start" << endl;

    Graph <int> G;
    int nv = 50;
    for(int i = 0; i < nv; i++){
        G.addVertex(i);
    }
    vector<int> pairs;
    int ne = 500;
    for(int i = 0; i < ne; i++){
        int r1 = rand() % nv;
        int r2 = rand() % nv;
        while(r2 == r1){
            r2 = rand() % nv;
        }
        if(!isPair(pairs, r1, r2)){
            int cost = rand() % 100;
            G.addEdge(r1, r2, 0, cost, 0);
            addPair(pairs, r1, r2);
        }
    }

    GraphPrintInfo(G);

    /*
    for(Vertex<int>* vertex : G.getVertexSet()){
        for(Edge<int>* edge : vertex->getOutgoing()){
            EdgePrintInfo(edge);
        }
    }
     */

    /*
    G.dijkstraShortestPath(0);

    vector <int> path1 = G.getPath(0, 1);


    for(int i = 0; i < path1.size(); i++){
        cout << path1[i] << endl;
    }


    for(Vertex<int>* vertex : G.getVertexSet()){
        VertexPrintInfo(vertex);
        cout << ": " << vertex->getDist() << endl;
    }

    G.dijkstraShortestPath(1);

    cout << "RIFT" << endl;


    for(Vertex<int>* vertex : G.getVertexSet()){
        VertexPrintInfo(vertex);
        cout << ": " << vertex->getDist() << endl;
    }
     */

    cout << G.getVertexSet().size() << endl;

    cout << "D DAY" << endl;

    vector<int> ip = {0, 27, 2, 9, 48, 4, 33};
    //Graph<int> mini_g = dijkstraInterestPoints(G, ip);
    Graph<int> mini_g = G.generateInterestPointsGraph(ip);
    GraphPrintInfo(mini_g);
    cout << mini_g.getVertexSet().size() << endl;

}