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

vector<double> countTimes(Graph<unsigned> ip_graph, unsigned bakery, vector<Van> vans) {
    vector<double> times;
    auto start = chrono::high_resolution_clock::now();; // Record start time
    auto finish = chrono::high_resolution_clock::now(); // Record end time
    start = chrono::high_resolution_clock::now();
    ip_graph.dividingClustersBrute(vans, bakery);
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
    map = Graph<unsigned>();
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
        //ip_map.viewGraphPath(path.first, path.first, true, true, true, true);
    }
}

void calculateFleetRoutes(unsigned bakery, Graph<unsigned> &main_map, Graph<unsigned> &ip_map) {
    Menu route_menu("Calcular rotas para frotas");
    route_menu.pushOption("Adicionar carrinha");
    route_menu.pushOption("Definir carrinhas predefinidas");
    route_menu.pushOption("Eliminar todas as carrinhas");
    route_menu.pushOption("Clusters + Greedy");
    route_menu.pushOption("Clusters + Brute Force");
    route_menu.pushOption("Voltar");
    unsigned choice, visit_time, capacity;
    double cost;
    vector<Van> vans;
    vector<unsigned> view_path;
    vector<unsigned> ips;
    std::vector<std::pair<Van, std::vector<Vertex<unsigned> *>>> van_paths;
    for (Vertex<unsigned> *v: ip_map.getVertexSet()) {
        ips.push_back(v->getInfo());
    }

    while (true) {
        choice = route_menu.chooseOption();
        switch (choice) {
            case 0:
                visit_time = getUnsigned("Tempo de visita");
                capacity = getUnsigned("Capacidade");
                vans.push_back(Van(vans.size(), visit_time, capacity));
                break;
            case 1:
                vans.clear();
                vans.push_back(Van(1, 2, 10));
                vans.push_back(Van(2, 3, 10));
                vans.push_back(Van(3, 5, 10));
                break;
            case 2:
                vans.clear();
                break;
            case 3:
                van_paths = ip_map.dividingClustersGreedy(vans, bakery);
                ip_map.followVansPath(van_paths, bakery);
                view_path = ip_map.vanPathtoViewPath(van_paths, bakery);
                ip_map.viewGraphPath(view_path, ips, true, true, true, true, true);
                //main_map.viewGraphPathIP(ip_map, view_path, true, true, true, true, true);
                break;
            case 4:
                van_paths = ip_map.dividingClustersBrute(vans, bakery);
                ip_map.followVansPath(van_paths, bakery);
                view_path = ip_map.vanPathtoViewPath(van_paths, bakery);
                ip_map.viewGraphPath(view_path, ips, true, true, true, true, true);
                //main_map.viewGraphPathIP(ip_map, view_path, true, false, true, true, true);//TODO: WORKING
                break;
            case 5:
                return;
            default:
                break;
        }

    }
}

void mainMenu() {
    Graph<unsigned> main_map = Graph<unsigned>();
    Graph<unsigned> ip_map = Graph<unsigned>();
    vector<unsigned> ids, ips;
    bool modif = false;
    unsigned bakery = UINF, client, vertex, time, hour, tolerance, max_tolerance;
    auto start = chrono::high_resolution_clock::now();
    auto finish = chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    Menu main_menu("Menu principal");
    main_menu.pushOption("Carregar mapa");
    main_menu.pushOption("Visualizar mapa");
    main_menu.pushOption("Definir padaria e propriedades");
    main_menu.pushOption("Definir cliente");
    main_menu.pushOption("Inserir impedimentos de transito");
    main_menu.pushOption("Calcular rotas");
    main_menu.pushOption("Calcular rotas para frotas");
    main_menu.pushOption("Analisar conectividade dos pontos de interesse");
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
                main_map.viewGraphIP(ip_map);
                break;
            case 2: // Definir padaria
                bakery = getUnsigned("Indice da padaria");
                while (main_map.findVertex(bakery) == nullptr) {
                    cout << "Indice inexistente! Insira um índice existente.";
                    bakery = getUnsigned("Indice da padaria");
                }
                time = getUnsigned("Velocidade media das carrinhas (km/h)");
                main_map.setVelocity(time * (1000 / 60));
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
                } else {
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
            case 6: // Calcular Rotas para Frotas
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
                calculateFleetRoutes(bakery, main_map, ip_map);
                break;
            case 7:
                ips = ids;
                ips.push_back(bakery);
                start = chrono::high_resolution_clock::now(); // Record start time
                if (main_map.checkConectivity(ips)) {
                    cout
                            << "Passou teste da conectividade.\nTodos os pontos de interesse fazem parte do mesmo componente fortemente conexo.\n";
                } else
                    cout
                            << "Não passou no teste da conectividade.\nOs pontos de interesse fazem parte do mesmo componente fortemente conexo.\n";
                finish = chrono::high_resolution_clock::now(); // Record end time
                elapsed = finish - start;
                break;
            case 8: // Reset para dados predefinidos
                ids.clear();
                ips.clear();
                main_map = Graph<unsigned>();
                ip_map = Graph<unsigned>();
                cout << "Todos os valores foram apagados.\n";
                cout << "A importar o mapa, pode demorar algum tempo...\n";
                main_map.importGraph("../resources/Porto/porto_strong_nodes_xy.txt",
                                     "../resources/Porto/porto_strong_edges.txt", false);
                cout << "A definir padaria e clientes...\n";
                ids = {174, 9, 11, 26, 26806, 26809, 26820, 47, 62};
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
            case 9: // Reset
                bakery = UINF;
                ids.clear();
                ips.clear();
                main_map = Graph<unsigned>();
                ip_map = Graph<unsigned>();
                modif = true;
                cout << "Todos os valores foram apagados.\n";
                break;
            case 10: // Sair
                quit = true;
                break;
            default:
                break;
        }
    }
}

int main() {
    mainMenu();
    return 0;

    // ------------- TIME TESTS ---------------
    Graph<unsigned> main_g;
    cout << "Importing graph...\n";
    main_g.importGraph("../resources/Porto/porto_strong_nodes_xy.txt", "../resources/Porto/porto_strong_edges.txt",
                       false);
    Graph<unsigned> ip_g;
    vector<double> bf_times;
    vector<double> times;
    vector<unsigned> ip_ids = {174, 9, 11, 26, 26806, 26809, 26820, 47, 16, 17, 21, 26, 27, 28, 30, 37, 39, 41, 47,
                               50, 53, 56, 57, 60, 91, 63, 66, 67, 69, 71, 73, 74, 75, 76, 80, 81, 82, 83, 84, 88, 89};
    vector<unsigned> input_ids;
    main_g.setEarlyTime(5);
    main_g.setStartTime(420);
    main_g.setVelocity(800);
    main_g.setVisitTime(5);

    vector<Van> vans;
    vans.push_back(Van(1, 2, 10));
    vans.push_back(Van(2, 3, 10));
    vans.push_back(Van(3, 5, 10));

    main_g.findVertex(9)->setTimes(430, 5, 10);
    main_g.findVertex(26)->setTimes(435, 5, 10);
    main_g.findVertex(26806)->setTimes(440, 5, 10);
    main_g.findVertex(26809)->setTimes(445, 5, 15);
    main_g.findVertex(26820)->setTimes(450, 5, 20);
    main_g.findVertex(47)->setTimes(455, 5, 10);
    main_g.findVertex(62)->setTimes(460, 5, 10);
    main_g.findVertex(11)->setTimes(470, 5, 8);
    main_g.findVertex(16)->setTimes(475, 5, 10);
    main_g.findVertex(17)->setTimes(480, 5, 10);
    main_g.findVertex(21)->setTimes(485, 5, 10);
    main_g.findVertex(26)->setTimes(490, 5, 10);
    main_g.findVertex(27)->setTimes(495, 5, 10);
    main_g.findVertex(28)->setTimes(500, 5, 10);
    main_g.findVertex(30)->setTimes(505, 5, 10);
    main_g.findVertex(37)->setTimes(510, 5, 10);
    main_g.findVertex(39)->setTimes(515, 5, 10);
    main_g.findVertex(41)->setTimes(520, 5, 10);
    main_g.findVertex(47)->setTimes(525, 5, 10);
    main_g.findVertex(50)->setTimes(530, 5, 10);
    main_g.findVertex(53)->setTimes(535, 5, 10);
    main_g.findVertex(56)->setTimes(540, 5, 10);
    main_g.findVertex(57)->setTimes(545, 5, 10);
    main_g.findVertex(60)->setTimes(550, 5, 10);
    main_g.findVertex(91)->setTimes(555, 5, 10);
    main_g.findVertex(63)->setTimes(560, 5, 10);
    main_g.findVertex(67)->setTimes(565, 5, 10);
    main_g.findVertex(69)->setTimes(570, 5, 10);
    main_g.findVertex(71)->setTimes(575, 5, 10);
    main_g.findVertex(73)->setTimes(580, 5, 10);
    main_g.findVertex(74)->setTimes(585, 5, 10);
    main_g.findVertex(75)->setTimes(590, 5, 10);
    main_g.findVertex(76)->setTimes(595, 5, 10);
    main_g.findVertex(80)->setTimes(600, 5, 10);
    main_g.findVertex(81)->setTimes(605, 5, 10);
    main_g.findVertex(82)->setTimes(610, 5, 10);
    main_g.findVertex(83)->setTimes(615, 5, 10);
    main_g.findVertex(84)->setTimes(620, 5, 10);
    main_g.findVertex(88)->setTimes(625, 5, 10);
    main_g.findVertex(89)->setTimes(630, 5, 10);

/*
    main_g.findVertex(9)->setTimes(430, 5, 300);
    main_g.findVertex(26)->setTimes(430, 5, 300);
    main_g.findVertex(26806)->setTimes(430, 5, 300);
    main_g.findVertex(26809)->setTimes(430, 5, 300);
    main_g.findVertex(26820)->setTimes(430, 5, 300);
    main_g.findVertex(47)->setTimes(430, 5, 300);
    main_g.findVertex(62)->setTimes(430, 5, 300);
    main_g.findVertex(11)->setTimes(430, 5, 300);
    main_g.findVertex(16)->setTimes(430, 5, 300);
    main_g.findVertex(17)->setTimes(430, 5, 300);
    main_g.findVertex(21)->setTimes(430, 5, 300);
    main_g.findVertex(26)->setTimes(430, 5, 300);
    main_g.findVertex(27)->setTimes(430, 5, 300);
    main_g.findVertex(28)->setTimes(430, 5, 300);
    main_g.findVertex(30)->setTimes(430, 5, 300);
    main_g.findVertex(37)->setTimes(430, 5, 300);
    main_g.findVertex(39)->setTimes(430, 5, 300);
    main_g.findVertex(41)->setTimes(430, 5, 300);
    main_g.findVertex(47)->setTimes(430, 5, 300);
    main_g.findVertex(50)->setTimes(430, 5, 300);
    main_g.findVertex(53)->setTimes(430, 5, 300);
    main_g.findVertex(56)->setTimes(430, 5, 300);
    main_g.findVertex(57)->setTimes(430, 5, 300);
    main_g.findVertex(60)->setTimes(430, 5, 300);
*/

    cout << "Calculating times...\n";

    unsigned counter = 1;
    string nn_out_file = "bf_times_out.txt";
    ofstream nn_times_file(nn_out_file);
    for (auto id = ip_ids.begin() + 2; id <= ip_ids.end(); id++) {
        input_ids = vector<unsigned>(ip_ids.begin(), id);
        ip_g = main_g.generateInterestPointsGraph(input_ids);
        cout << counter << " clientes\n";
        times = countTimes(ip_g, 174, vans);
        nn_times_file << times.at(0) << endl;
        bf_times.push_back(times.at(0));
        counter++;
    }


    /*
    string nn_out_file = "bf_times_out.txt";
    cout << "Writing results to " << nn_out_file << endl;
    ofstream nn_times_file(nn_out_file);
    for (auto time : bf_times) {
        nn_times_file << time << endl;
    }
*/
    nn_times_file.close();
    cout << "Done!\n";

    return 0;
    // ------------- END TIME TESTS ---------------
}
