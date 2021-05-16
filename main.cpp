#include <iostream>
#include <string>
#include "Graph.h"
#include <math.h>
#include <ctime>

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
template<class T>
Graph<T> importGraph(string vertex_filename, string edges_filename) {
    ifstream vertex_file(vertex_filename);
    ifstream edge_file(edges_filename);
    Graph<unsigned> graph;
    unsigned n_vertex;
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
        graph.addVertex(id, x, y);
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
        graph.addEdge(orig, dest);
    }
    return graph;
}

void viewGraph(Graph<unsigned> graph) {
    // Instantiate GraphViewer
    GraphViewer gv;

    // Set coordinates of window center
    gv.setCenter(sf::Vector2f(300, 300));

    for (auto vertex : graph.getVertexSet()) {
        gv.addNode(vertex->getInfo(), sf::Vector2f(vertex->getX(), vertex->getY()));
    }

    unsigned counter = 0;
    for (auto vertex : graph.getVertexSet()) {
        for (auto edge : vertex->getOutgoing()) {
            gv.addEdge(counter, gv.getNode(vertex->getInfo()), gv.getNode(edge->getDest()->getInfo()), GraphViewer::Edge::DIRECTED);
            counter++;
        }
    }

    // Create window
    gv.createWindow(600, 600);

    // Join viewer thread (blocks till window closed)
    gv.join();
}
 */

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

    /*
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

    /*
    cout << G.getVertexSet().size() << endl;

    cout << "D DAY" << endl;

    vector<int> ip = {0, 27, 2, 9, 48, 4, 33};
    //Graph<int> mini_g = dijkstraInterestPoints(G, ip);
    Graph<int> mini_g = G.generateInterestPointsGraph(ip);
    GraphPrintInfo(mini_g);
    cout << mini_g.getVertexSet().size() << endl;
     */

    Graph<unsigned> g;
    Graph<unsigned> gg;
    cout << "Importing graph..." << endl;
    time_t start_time = time(NULL);
    g.importGraph("../resources/Porto/porto_strong_nodes_xy.txt", "../resources/Porto/porto_strong_edges.txt");
    cout << "Calculating scc..." << endl;
    cout << "STRONG Detected " << g.dfsConnectivity().size() << " scc's" << endl;
    cout << "Importing graph..." << endl;
    gg.importGraph("../resources/Porto/porto_full_nodes_xy.txt", "../resources/Porto/porto_full_edges.txt");
    cout << "Calculating scc..." << endl;
    cout << "FULL Detected " << gg.dfsConnectivity().size() << " scc's" << endl;
    //g.importGraph("../resources/Espinho/espinho_strong_nodes_xy.txt", "../resources/Espinho/espinho_strong_edges.txt");
    time_t end_time = time(NULL);
    //cout << "Finished importing graph in " << end_time - start_time << " s" << endl;
    //vector<unsigned> ids {53619, 21329, 48150, 26850};
    //vector<unsigned> ids {8932, 13373};
    //cout << "Running Dijkstra..." << endl;
    //Graph<unsigned> minig = g.generateInterestPointsGraph(ids);
    //GraphPrintInfo(g);

    //minig.viewGraph();
    //g.viewGraphPath(minig);
}
