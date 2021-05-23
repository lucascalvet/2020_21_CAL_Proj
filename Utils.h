#ifndef INC_2020_21_CAL_PROJ_UTILS_H
#define INC_2020_21_CAL_PROJ_UTILS_H

#include <string>

double calculateDist(double x0, double y0, double x1, double y1);

double calculateDistHaversine(double lat0, double long0, double lat1, double long1);

void stopConsole();

unsigned getUnsigned(const std::string &question);

#endif //INC_2020_21_CAL_PROJ_UTILS_H
