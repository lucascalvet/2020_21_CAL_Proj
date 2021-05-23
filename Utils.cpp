#include "Utils.h"
#include <cmath>
#include <string>
#include <iostream>
#include <limits>

using namespace std;
#define EARTH_R 6371.0

double calculateDist(double x0, double y0, double x1, double y1) {
    return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
}

double calculateDistHaversine(double long0, double lat0, double long1, double lat1) {
    double rad_lat0 = lat0 * (M_PI / 180.0);
    double rad_lat1 = lat1 * (M_PI / 180.0);
    double rad_lat_diff = (lat1 - lat0) * (M_PI / 180.0);
    double rad_long_diff = (long1 - long0) * (M_PI / 180.0);

    double trig = pow(sin(rad_lat_diff / 2), 2) +
               pow(sin(rad_long_diff / 2), 2) *
               cos(rad_lat0) * cos(rad_lat1);
    double dist_km = EARTH_R * 2 * asin(sqrt(trig));
    return dist_km * 1000.0;
}

void stopConsole() {
    cout << endl; //formatting console
    cout << "Press enter to continue...";
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

unsigned getUnsigned(const string& question) {
    unsigned choice;
    cout << question << ": ";
    cin >> choice;
    while (cin.fail() || cin.peek() != '\n') {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Valor inválido!\n";
        cout << question << ": ";
        cin >> choice;
    }
    return choice;
}
