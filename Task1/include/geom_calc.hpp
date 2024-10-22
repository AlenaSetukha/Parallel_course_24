#ifndef _GEOM_CALC_H_
#define _GEOM_CALC_H_

#include <iostream>
#include <vector>
#include "Point.hpp"

/**
 * Функция, вычисляющая площадь выпуклого
 * N-угольника по координатам вершин(x, y).
 */

double get_Area(const std::vector<Point>& vertices);

double scal_prod(const std::vector<std::vector<double>>& a, 
                 const std::vector<std::vector<double>>& b, 
                 const double h1, const double h2);
                 
double get_normC(std::vector<std::vector<double>> a, 
                 std::vector<std::vector<double>> b);

#endif