#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>
#include "Graph.h"
#include "Menu.h"
#include "Utils.h"

using namespace std;

template<class T>
void VertexPrintInfo(Vertex<T> *v) { cout << "V" << v->getInfo(); }

template<class T>
void EdgePrintInfo(Edge<T> *e) {
    VertexPrintInfo(e->getOrigin());
    cout << "--" << e->getCost() << "-->";
    VertexPrintInfo(e->getDest());
    cout << endl;
}

template<class T>
void GraphPrintInfo(Graph<T> g) {
    for (Vertex<T> *vertex : g.getVertexSet()) {
        for (Edge<T> *edge : vertex->getOutgoing()) {
            EdgePrintInfo(edge);
        }
        /*
        for(Edge<T>* edge : vertex->getIncoming()){
            EdgePrintInfo(edge);
        }
         */
    }
}

template<class T>
void PrintVector(vector<T> vec, string title) {
    cout << title << ": ";
    for (int i = 0; i < vec.size(); i++) {
        cout << vec[i];
        if (i < vec.size() - 1) {
            cout << ", ";
        } else cout << endl;
    }
}

vector<double> countTimes(Graph<unsigned> ip_graph , unsigned bakery) {
    vector<double> times;

    auto start = chrono::high_resolution_clock::now();; // Record start time
    auto finish = chrono::high_resolution_clock::now(); // Record end time
    start = chrono::high_resolution_clock::now();
    ip_graph.heldKarp(bakery);
    finish = chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    times.push_back(elapsed.count());
    return times;
}

void importMap(Graph<unsigned> &map) {
    Menu map_menu("Carregar mapa");
    map_menu.pushOption("Porto");
    map_menu.pushOption("Espinho");
    map_menu.pushOption("Importar de ficheiro (x, y)");
    map_menu.pushOption("Importar de ficheiro (lat, long)");
    map_menu.pushOption("Voltar");

    string nodes_filename, edges_filename;
    bool lat_lng = false;
    unsigned choice;
    choice = map_menu.chooseOption();
    switch (choice) {
        case 0: // Porto
            nodes_filename = "../resources/Porto/porto_strong_nodes_xy.txt";
            edges_filename = "../resources/Porto/porto_strong_edges.txt";
            break;
        case 1: // Espinho
            nodes_filename = "../resources/Espinho/espinho_strong_nodes_xy.txt";
            edges_filename = "../resources/Espinho/espinho_strong_edges.txt";
            break;
        case 2: // Importar de ficheiro (x, y)
        case 3: // Importar de ficheiro (lat, long)
            cout << "Ficheiro de vertices: ";
            cin >> nodes_filename;
            while (cin.fail() || cin.peek() != '\n') {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Input invalido! Insira um caminho de ficheiro valido.\n";
                cout << "Caminho do ficheiro de vertices: ";
                cin >> nodes_filename;
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Ficheiro de arestas: ";
            cin >> edges_filename;
            while (cin.fail() || cin.peek() != '\n') {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Input invalido! Insira um caminho de ficheiro valido.\n";
                cout << "Caminho do ficheiro de arestas: ";
                cin >> edges_filename;
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            break;
        case 4: // Voltar
            return;
        default:
            break;
    }
    if (choice == 3) lat_lng = true;
    cout << "A importar mapa, pode demorar algum tempo...\n";
    map.importGraph(nodes_filename, edges_filename, lat_lng);
    cout << "Mapa importado com sucesso.\n";
}

void calculateRoutes(unsigned bakery, Graph<unsigned> &main_map, Graph<unsigned> &ip_map) {
    Menu route_menu("Calcular rotas");
    route_menu.pushOption("Held-Karp (sem contagem de horas)");
    route_menu.pushOption("Nearest Neighbour (sem contagem de horas)");
    route_menu.pushOption("Nearest Neighbour (com contagem de horas)");
    route_menu.pushOption("Brute Force (com contagem de horas)");
    route_menu.pushOption("Voltar");
    unsigned choice;
    double cost;
    while (true) {
        choice = route_menu.chooseOption();
        if (choice != 4) {
            cout << "A correr ";
        }
        auto start = chrono::high_resolution_clock::now(); // Record start time
        switch (choice) {
            case 0: // Held-Karp (sem contagem de horas)
                cout << "Held-Karp (sem contagem de horas)...\n";
                ip_map.heldKarp(bakery);
                break;
            case 1: // Nearest Neighbour (sem contagem de horas)
                cout << "Nearest Neighbour (sem contagem de horas)...\n";
                ip_map.nearestNeighbour(bakery);
                break;
            case 2: // Nearest Neighbour (com contagem de horas)
                cout << "Nearest Neighbour (com contagem de horas)...\n";
                ip_map.nearestNeighbourTimes(bakery);
                break;
            case 3: // Brute Force (com contagem de horas)
                cout << "Brute Force (com contagem de horas)...\n";
                cost = ip_map.bruteForceTimes(bakery);
                break;
            case 4: // Voltar
            default:
                return;
        }
        auto finish = chrono::high_resolution_clock::now(); // Record end time
        std::chrono::duration<double> elapsed = finish - start;
        if (choice > 1) {
            ip_map.printTimes();
            if (choice == 3) {
                cout << "Função de custo: " << cost << '\n';
            }
        }
        cout << "Resultado obtido em " << elapsed.count() << "s.\n";
        pair<vector<unsigned>, double> path = ip_map.getPath(bakery, bakery);
        cout << "Distância percorrida: " << path.second << '\n';
        main_map.viewGraphPathIP(ip_map, path.first, true, false, true, true);
        ip_map.viewGraphPath(path.first, path.first, true, true, true, true);
    }
}

void mainMenu() {
    Graph<unsigned> main_map = Graph<unsigned>();
    Graph<unsigned> ip_map = Graph<unsigned>();
    vector<unsigned> ids, ips;
    bool modif = false;
    unsigned bakery = UINF, client, vertex, time, hour, tolerance, max_tolerance;
    Menu main_menu("Menu principal");
    main_menu.pushOption("Carregar mapa");
    main_menu.pushOption("Visualizar mapa");
    main_menu.pushOption("Definir padaria e propriedades");
    main_menu.pushOption("Definir cliente");
    main_menu.pushOption("Inserir impedimentos de transito");
    main_menu.pushOption("Calcular rotas");
    main_menu.pushOption("Reset para dados predefinidos");
    main_menu.pushOption("Reset");
    main_menu.pushOption("Sair");
    unsigned choice;
    bool quit = false;
    while (!quit) {
        choice = main_menu.chooseOption();
        switch (choice) {
            case 0: // Carregar mapa
                importMap(main_map);
                break;
            case 1: // Visualizar mapa
                ips = ids;
                if (bakery != UINF) ips.push_back(bakery);
                if (modif) {
                    cout
                            << "A correr o algoritmo de Dijkstra para calcular os caminhos mais curtos entre os pontos de interesse...";
                    ip_map = main_map.generateInterestPointsGraph(ips);
                    cout << "Concluido.\n";
                    modif = false;
                }
                cout << "A abrir GraphViewer...\n";
                main_map.viewGraph();
                break;
            case 2: // Definir padaria
                bakery = getUnsigned("Indice da padaria");
                while (main_map.findVertex(bakery) == nullptr) {
                    cout << "Indice inexistente! Insira um índice existente.";
                    bakery = getUnsigned("Indice da padaria");
                }
                time = getUnsigned("Velocidade media das carrinhas (km/h)");
                main_map.setVelocity(time * (1000/60));
                time = getUnsigned("Hora de saida da padaria (min desde as 0:00)");
                main_map.setEarlyTime(time);
                time = getUnsigned("Tempo permitido de antecedencia (min)");
                main_map.setEarlyTime(time);
                time = getUnsigned("Tempo de visita (min)");
                main_map.setVisitTime(time);
                modif = true;
                break;
            case 3: // Definir cliente
                client = getUnsigned("Indice do cliente");
                while (main_map.findVertex(client) == nullptr) {
                    cout << "Indice inexistente! Insira um indice existente.\n";
                    client = getUnsigned("Indice do cliente");
                }
                if (client == bakery) {
                    cout << "A padaria já foi definida no vértice " << bakery << "! Operação abortada.\n";
                    stopConsole();
                    break;
                }
                hour = getUnsigned("Hora pretendida de entrega (min desde as 0:00)");
                tolerance = getUnsigned("Tolerancia de tempo (min)");
                max_tolerance = getUnsigned("Tolerancia maxima de tempo (min)");
                main_map.findVertex(client)->setTimes(hour, tolerance, max_tolerance);
                if (find(ids.begin(), ids.end(), client) == ids.end()) {
                    ids.push_back(client);
                    cout << "Cliente foi definido.\n";
                }
                else {
                    cout << "Cliente foi redefinido.\n";
                }
                break;
            case 4: // Inserir impedimentos de trânsito
                Vertex<unsigned> *v1;
                Vertex<unsigned> *v2;
                Edge<unsigned> *e;
                vertex = getUnsigned("Indice do vértice de saída");
                v1 = main_map.findVertex(vertex);
                if (v1 == nullptr) {
                    cout << "Índice inexistente!\n";
                    break;
                }
                vertex = getUnsigned("Indice do vértice de chegada");
                v2 = main_map.findVertex(vertex);
                if (v2 == nullptr) {
                    cout << "Índice inexistente!\n";
                    break;
                }
                e = main_map.findEdge(v1, v2);
                if (e == nullptr) {
                    cout << "Não existe qualquer estrada a ligar " << v1->getInfo() << " a " << v2->getInfo() << "!\n";
                    break;
                }
                v1->removeEdge(e);
                break;
            case 5: // Calcular rotas
                if (bakery == UINF) {
                    cout << "A Padaria nao esta definida!\n";
                    stopConsole();
                    break;
                }
                if (ids.empty()) {
                    cout << "Nenhum cliente definido!\n";
                    stopConsole();
                    break;
                }
                ips = ids;
                ips.push_back(bakery);
                if (modif) {
                    cout
                            << "A correr o algoritmo de Dijkstra para calcular os caminhos mais curtos entre os pontos de interesse...";
                    ip_map = main_map.generateInterestPointsGraph(ips);
                    cout << "Concluido.\n";
                    modif = false;
                }
                calculateRoutes(bakery, main_map, ip_map);
                break;
            case 6: // Reset para dados predefinidos
                ids.clear();
                ips.clear();
                main_map = Graph<unsigned>();
                ip_map = Graph<unsigned>();
                cout << "Todos os valores foram apagados.\n";
                cout << "A importar o mapa...\n";
                main_map.importGraph("../resources/Porto/porto_strong_nodes_xy.txt",
                                     "../resources/Porto/porto_strong_edges.txt", false);
                cout << "A definir padaria e clientes...\n";
                ids = {9, 11, 26, 26806, 26809, 26820, 47, 62};
                bakery = 174;
                main_map.setEarlyTime(5);
                main_map.setStartTime(420);
                main_map.setVelocity(800);
                main_map.setVisitTime(5);
                main_map.findVertex(9)->setTimes(430, 5, 10);
                main_map.findVertex(26)->setTimes(435, 5, 10);
                main_map.findVertex(26806)->setTimes(440, 5, 10);
                main_map.findVertex(26809)->setTimes(445, 5, 15);
                main_map.findVertex(26820)->setTimes(450, 5, 20);
                main_map.findVertex(47)->setTimes(455, 5, 10);
                main_map.findVertex(62)->setTimes(460, 5, 10);
                main_map.findVertex(11)->setTimes(470, 5, 8);
                modif = true;
                cout << "Concluido.\n";
                break;
            case 7: // Reset
                bakery = UINF;
                ids.clear();
                ips.clear();
                main_map = Graph<unsigned>();
                ip_map = Graph<unsigned>();
                modif = true;
                cout << "Todos os valores foram apagados.\n";
                break;
            case 8: // Sair
                quit = true;
                break;
            default:
                break;
        }
    }
}

int main() {
    //mainMenu();
    //return 0;

    // ------------- TIME TESTS ---------------
    Graph<unsigned> main_g;
    cout << "Importing graph...\n";
    main_g.importGraph("../resources/Porto/porto_strong_nodes_xy.txt", "../resources/Porto/porto_strong_edges.txt", false);
    Graph<unsigned> ip_g;
    vector<double> hk_times;
    vector<double> times;
    vector<unsigned> ip_ids = {174, 9, 11, 26, 26806, 26809, 26820, 47, 62};
    vector<unsigned> input_ids;
    main_g.setEarlyTime(5);
    main_g.setStartTime(420);
    main_g.setVelocity(800);
    main_g.setVisitTime(5);
    cout << "Calculating times...\n";
    unsigned counter = 1;
    for (auto id = ip_ids.begin() + 2; id <= ip_ids.end(); id++) {
        input_ids = vector<unsigned>(ip_ids.begin(), id);
        ip_g = main_g.generateInterestPointsGraph(input_ids);
        cout << counter << " clientes\n";
        times = countTimes(ip_g, 174);
        hk_times.push_back(times.at(0));
        counter++;
    }

    string out_file = "hk_times_out.txt";
    cout << "Writing results to " << out_file << endl;
    ofstream hk_times_file(out_file);
    for (auto time : hk_times) {
        hk_times_file << time << endl;
    }

    hk_times_file.close();
    cout << "Done!\n";

    return 0;
    // ------------- END TIME TESTS ---------------

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
    //g.importGraph("../resources/Porto/porto_strong_nodes_xy.txt", "../resources/Porto/porto_strong_edges.txt", false);
    g.importGraph("../resources/Porto/porto_strong_nodes_latlng.txt", "../resources/Porto/porto_strong_edges.txt", true);
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
    g.findVertex(9)->setTimes(100, 5, 10);
    g.findVertex(26)->setTimes(150, 5, 10);
    g.findVertex(26806)->setTimes(160, 5, 10);
    g.findVertex(26809)->setTimes(170, 5, 10);
    g.findVertex(26820)->setTimes(490, 5, 10);
    g.findVertex(47)->setTimes(500, 5, 10);
    g.findVertex(62)->setTimes(600, 5, 10);
    vector<unsigned> ids = {9, 11, 26, 26806, 26809, 26820, 47, 62};
    //vector<unsigned> ids {8932, 13373};
    cout << "Running Dijkstra..." << endl;
    Graph<unsigned> minig = g.generateInterestPointsGraph(ids);
    //Graph<unsigned> minig = gg.generateInterestPointsGraph(ids);
    //GraphPrintInfo(g);

    /*
    std::vector<Cluster<unsigned>> vc = minig.getClusters(11);

    for(Cluster<unsigned> c : vc){
        c.PrintInfo();
    }

    std::vector<Van> vv;
    vv.push_back(Van(1, 5, 10));
    //vv.push_back(Van(2, 40, 10));
    vv.push_back(Van(3, 2, 10));
    //vv.push_back(Van(4, 12, 10));
    //vv.push_back(Van(5, 1, 10));

    cout << "Dividing Clusters Greedy..." << endl;

    //vector<pair<Van, vector<Vertex<unsigned > *>>> van_pairs = minig.dividingClustersGreedy(vv, 11);
    vector<pair<Van, vector<Vertex<unsigned > *>>> van_pairs = minig.dividingClustersBrute(vv, 11);
    for(int i = 0; i < van_pairs.size(); i++){
        cout << "Van " << van_pairs[i].first.getID() << ": ";
        for(int j = 0; j < van_pairs[i].second.size(); j++){
            cout << van_pairs[i].second[j]->getInfo() << " ";
        }
        cout << "END" << endl;
    }

    minig.followVansPath(van_pairs, 11);

    return 0;
     */

    //GraphPrintInfo(minig);
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

    //vector<unsigned> ov = minig.getOverlapClients(26806);
    //vector<unsigned> ovt = minig.getOverlapClientsTravelling(26806);

    cout << "Running Nearest Neighbour Times..." << endl;
    minig.nearestNeighbourTimes(11);
    cout << "Getting path..." << endl;
    pair<vector<unsigned>, double> nnt_path = minig.getPath(11, 11);

    minig.printTimes();
    cout << "GraphViewer..." << endl;
    minig.viewGraphPath(nnt_path.first, ids, true, true, true, true, true);


    cout << "Running Brute Force Times..." << endl;
    if (minig.bruteForceTimes(11) == INF) cout << "Impossible path!\n";
    else {
        cout << "Getting path..." << endl;
        pair<vector<unsigned>, double> bf_path = minig.getPath(11, 11);

        minig.printTimes();

        cout << "GraphViewer..." << endl;
        minig.viewGraphPath(bf_path.first, ids, true, true, true, true, true);
    }

    return 0;

    cout << "Running Held-Karp..." << endl;
    minig.heldKarp(11);

    cout << "Getting path..." << endl;
    pair<vector<unsigned>, double> hk_path = minig.getPath(11, 11);

    cout << "GraphViewer..." << endl;

    /*
    if(minig.getVertexSet().size() != 0){
        g.viewGraphIP(minig);
    }
    */

    //PrintVector(ov, "OV[26806]");
    //PrintVector(ovt, "OVT[26806]");


    //minig.viewGraphPath(hk_path.first, ids, false, true);
    //g.viewGraph();

    //g.viewGraphPathIP(minig, hk_path.first, true, true);

    minig.viewGraphPath(hk_path.first, ids, true, true, true, true);

    cout << "Running Nearest Neighbour..." << endl;
    minig.nearestNeighbour(11);
    cout << "Getting path..." << endl;
    pair<vector<unsigned>, double> nn_path = minig.getPath(11, 11);

    cout << "GraphViewer..." << endl;
    minig.viewGraphPath(nn_path.first, ids, true, true, true, true);

    cout << "HK: " << hk_path.second << " vs NN: " << nn_path.second << endl;
    minig.printTimes();
    //g.viewGraphPathIP(minig, nn_path.first, true);


    double dist = calculateDistHaversine(-8.577122, 41.172792, -8.616987, 41.150174);
    //double dist = calculateDistHaversine(41.172792, -8.577122, 41.150174, -8.616987);
    double old_dist = calculateDist(6218.020297963754, -2100.583359242417, 1785.2345472782617, 414.4234916046262);
    cout << "XY DIST: " << old_dist << " vs LATLNG DIST: " << dist << endl;
}
