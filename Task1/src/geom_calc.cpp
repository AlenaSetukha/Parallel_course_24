#include <iostream>
#include <vector>
#include <numeric>
#include "Point.hpp"

double get_Area(const std::vector<Point>& vertices)
{
    double area = 0.0;
    int N = vertices.size();

    for (int i = 0; i < N; i++) {
        int i_next = (i + 1) % N;
        area += vertices[i].x * vertices[i_next].y;
        area -= vertices[i_next].x * vertices[i].y;
    }
    return std::abs(area) / 2.0;
}


double scal_prod(const std::vector<std::vector<double>>& a, 
                 const std::vector<std::vector<double>>& b, 
                 const double h1, const double h2) 
{
    double res = 0.0;
    int Ns = a.size(), Ms = a[0].size();
    
    for (int j = 1; j < Ns - 1; j++) { // y
        res += std::inner_product(a[j].begin() + 1, a[j].end() - 1, 
                                   b[j].begin() + 1, 0.0);
    }
    return res * h1 * h2;
}


double get_normC(std::vector<std::vector<double>> a, 
                 std::vector<std::vector<double>> b)
{
    if (a.size() != b.size() || a.empty() || b.empty() || a[0].size() != b[0].size()) {
        throw std::invalid_argument("Matrices size: get_normC");
    }
    double res = 0.0;
    int Ns = a.size(), Ms = a[0].size();
    for (int j = 1; j < Ns - 1; j++) { // y
        for (int i = 1; i < Ms - 1; i++) { // x
            double diff = std::abs(a[j][i] - b[j][i]);
            res = std::max(res, diff);
        }
    }
    return res;
}
