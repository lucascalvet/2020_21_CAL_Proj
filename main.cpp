#include <iostream>
#include "Graph.h"
#include <ctime>

using namespace std;

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

template <class T>
void PrintVector(vector<T> vec, string title) {
    cout << title << ": ";
    for(int i = 0; i < vec.size(); i++){
        cout << vec[i];
        if(i < vec.size() - 1){
            cout << ", ";
        }
        else cout << endl;
    }
}

int main(){
    cout << "Start" << endl;



    Graph<unsigned> g;
    Graph<unsigned> gg;
    /*
    cout << "Importing graph..." << endl;
    time_t start_time = time(NULL);
    g.importGraph("../resources/Porto/porto_strong_nodes_xy.txt", "../resources/Porto/porto_strong_edges.txt");
    cout << "Calculating scc..." << endl;
    cout << "STRONG Detected " << g.dfsConnectivity().size() << " scc's" << endl;
     */


    cout << "Importing graph..." << endl;
    g.importGraph("../resources/Porto/porto_strong_nodes_xy.txt", "../resources/Porto/porto_strong_edges.txt", false);
    //gg.importGraph("../resources/Porto/porto_full_nodes_xy.txt", "../resources/Porto/porto_full_edges.txt");
    //cout << "Calculating scc..." << endl;
    //cout << "FULL Detected " << gg.dfsConnectivity().size() << " scc's" << endl;
    //g.importGraph("../resources/Espinho/espinho_strong_nodes_xy.txt", "../resources/Espinho/espinho_strong_edges.txt");
    //time_t end_time = time(NULL);
    //cout << "Finished importing graph in " << end_time - start_time << " s" << endl;
    g.setEarlyTime(5);
    g.setStartTime(0);
    g.setVelocity(100);
    g.setVisitTime(5);
    g.findVertex(9)->setTimes(400, 5, 10);
    g.findVertex(26)->setTimes(100, 5, 10);
    g.findVertex(26806)->setTimes(200, 5, 10);
    g.findVertex(26809)->setTimes(250, 5, 10);
    g.findVertex(26820)->setTimes(500, 5, 10);
    g.findVertex(47)->setTimes(600, 5, 10);
    g.findVertex(62)->setTimes(650, 5, 10);
    vector<unsigned> ids = {9, 11, 26, 26806, 26809, 26820, 47, 62};
    //vector<unsigned> ids {8932, 13373};
    cout << "Running Dijkstra..." << endl;
    Graph<unsigned> minig = g.generateInterestPointsGraph(ids);
    //Graph<unsigned> minig = gg.generateInterestPointsGraph(ids);
    //GraphPrintInfo(g);

    GraphPrintInfo(minig);
    //minig.viewGraph();

    /*
    minig.findVertex(9)->setTimes(20, 5, 10);
    minig.findVertex(26)->setTimes(40, 5, 10);
    minig.findVertex(26806)->setTimes(50, 5, 10);
    minig.findVertex(26809)->setTimes(60, 5, 10);
    minig.findVertex(26820)->setTimes(110, 5, 10);
    minig.findVertex(47)->setTimes(125, 5, 10);
    minig.findVertex(62)->setTimes(140, 5, 10);
     */

    vector<unsigned> ov = minig.getOverlapClients(26806);
    vector<unsigned> ovt = minig.getOverlapClientsTravelling(26806);

    cout << "Running Held-Karp..." << endl;
    minig.heldKarp(11);

    cout << "Getting path..." << endl;
    pair<vector<unsigned> , double> hk_path = minig.getPath(11, 11);



    cout << "GraphViewer..." << endl;

    /*
    if(minig.getVertexSet().size() != 0){
        g.viewGraphIP(minig);
    }
    */

    PrintVector(ov, "OV[26806]");
    PrintVector(ovt, "OVT[26806]");


    //minig.viewGraphPath(hk_path.first, ids, false, true);
    //g.viewGraph();

    //g.viewGraphPathIP(minig, hk_path.first, true, true);

    minig.viewGraphPath(hk_path.first, ids, true, true, true);

    cout << "Running Nearest Neighbour..." << endl;
    //minig.nearestNeighbour(11);
    minig.nearestNeighbourTimes(11);
    cout << "Getting path..." << endl;
    pair<vector<unsigned> , double> nn_path = minig.getPath(11, 11);

    cout << "GraphViewer..." << endl;
    minig.viewGraphPath(nn_path.first, ids, true, true, true);

    cout << "HK: " << hk_path.second << " vs NN: " << nn_path.second << endl;
    minig.printTimes();
    //g.viewGraphPathIP(minig, nn_path.first, true);


    double dist = calculateDistHaversine(-8.577122, 41.172792, -8.616987, 41.150174);
    //double dist = calculateDistHaversine(41.172792, -8.577122, 41.150174, -8.616987);
    double old_dist = calculateDist(6218.020297963754,-2100.583359242417, 1785.2345472782617,414.4234916046262);
    cout << "XY DIST: " << old_dist << " vs LATLNG DIST: " << dist << endl;
}
