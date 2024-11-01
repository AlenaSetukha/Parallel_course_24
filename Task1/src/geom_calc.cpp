#include <iostream>
#include <vector>
#include <numeric>
#include <omp.h>
#include <cmath>
 

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
    if (a.empty() || b.empty() || a.size() != b.size() || a[0].empty()) {
        throw std::invalid_argument("Invalid input vectors");
    }
    double res = 0.0;
    int Ns = a.size(), Ms = a[0].size();
    
    for (int j = 1; j < Ns - 1; j++) { // y
        res += std::inner_product(a[j].begin() + 1, a[j].end() - 1, 
                                  b[j].begin() + 1, 0.0);
    }
    return res * h1 * h2;
}



double scal_prod_OMP(const std::vector<std::vector<double>>& a, 
                 const std::vector<std::vector<double>>& b, 
                 const int num_threads,
                 const double h1, const double h2) 
{
    if (a.empty() || b.empty() || a.size() != b.size() || a[0].empty()) {
        throw std::invalid_argument("Invalid input vectors");
    }
    double res = 0.0;
    int Ns = a.size(), Ms = a[0].size();
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for reduction(+ : res)
    for (int j = 1; j < Ns - 1; j++) { // y
        res += std::inner_product(a[j].begin() + 1, a[j].end() - 1, 
                                  b[j].begin() + 1, 0.0);
    }
    return res * h1 * h2;
}













double get_normC(const std::vector<std::vector<double>> a, 
                 const std::vector<std::vector<double>> b)
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


double get_normC_OMP(const std::vector<std::vector<double>> a, 
                const std::vector<std::vector<double>> b,
                const int num_threads)
{
    if (a.size() != b.size() || a.empty() || b.empty() || a[0].size() != b[0].size()) {
        throw std::invalid_argument("Matrices size: get_normC");
    }
    double res = 0.0;
    int Ns = a.size(), Ms = a[0].size();
    #pragma omp parallel
    {
        double local_res = 0.0;
        #pragma omp for
        for (int j = 1; j < Ns - 1; j++) { // y
            for (int i = 1; i < Ms - 1; i++) { // x
                double diff = std::abs(a[j][i] - b[j][i]);
                if (diff > local_res) {
                    local_res = diff;
                }
            }
        }
        // Объединение локальных результатов
        #pragma omp critical
        {
            if (local_res > res) {
                res = local_res;
            }
        }
    }
    return res;
}
