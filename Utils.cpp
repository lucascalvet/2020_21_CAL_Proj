#include "Utils.h"
#include <cmath>

using namespace std;
#define EARTH_R pow(6371.0, 3.0)

double calculateDist(double x0, double y0, double x1, double y1) {
    return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
}

double calculateDistHaversine(double long0, double lat0, double long1, double lat1) {
    double rad_lat0 = lat0 * (M_PI / 180.0);
    double rad_lat1 = lat1 * (M_PI / 180.0);
    double rad_lat_diff = (lat1 - lat0) * (M_PI / 180.0);
    double rad_long_diff = (long1 - long0) * (M_PI / 180.0);
    /*
    double trig = sin(rad_lat_diff / 2.0) * sin(rad_lat_diff / 2.0) +
                  cos(rad_lat0) * cos(rad_lat1) * sin(rad_long_diff / 2.0) * sin(rad_long_diff / 2.0);
    return EARTH_R * 2 * atan2(sqrt(trig), sqrt(1 - trig));
     */
    double trig = pow(sin(rad_lat_diff / 2), 2) +
               pow(sin(rad_long_diff / 2), 2) *
               cos(rad_lat0) * cos(rad_lat1);
    double dist_km = 6371.0 * 2 * asin(sqrt(trig));
    return dist_km * 1000.0;
}
