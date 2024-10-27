#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "Point.hpp"
#include "coeff_calc.hpp"
#include "geom_calc.hpp"

static inline double sqr(const double x)
{
    return x * x;
}

void get_matrixA(const std::vector<std::vector<double>>& a,
        const std::vector<std::vector<double>>& b,
        const std::vector<std::vector<double>>& w,
        const double h1, const double h2,
        std::vector<std::vector<double>>& matrix)
{
    double diff1, diff2;
    int Ns = w.size(), Ms = w[0].size();
    for (int j = 1; j < Ns - 1; j++)//y
    {
        for (int i = 1; i < Ms - 1; i++)//x
        {
            diff1 = a[j][i + 1] * (w[j][i + 1] - w[j][i]) -
                        a[j][i] * (w[j][i] - w[j][i - 1]);
            diff2 = b[j + 1][i] * (w[j + 1][i] - w[j][i]) -
                        b[j][i] * (w[j][i] - w[j - 1][i]);
            matrix[j][i] = -diff1 / sqr(h1) - diff2 / sqr(h2);
        }
    }
}



void get_matrixB(const std::vector<std::vector<Point>>& grid,
        std::vector<std::vector<double>>& matrix)
{
    int Ns = grid.size(), Ms = grid[0].size();
    double h1 = std::abs(grid[0][0].x - grid[0][1].x);
    double h2 = std::abs(grid[0][0].y - grid[1][0].y);
    for (int j = 1; j < Ns - 1; j++)//y
    {
        for (int i = 1; i < Ms - 1; i++)//x
        {
            matrix[j][i] = get_fCoeff(grid[j][i].x - 0.5 * h1, grid[j][i].x + 0.5 * h1,
                                      grid[j][i].y + 0.5 * h2, grid[j][i].y - 0.5 * h2);
        }
    }
}




void get_aCoeffMatrix(const std::vector<std::vector<Point>>& grid,
        const double eps,
        std::vector<std::vector<double>>& matrix)
{
    int Ns = grid.size(), Ms = grid[0].size();
    double h1 = std::abs(grid[0][0].x - grid[0][1].x);
    double h2 = std::abs(grid[0][0].y - grid[1][0].y);
    for (int j = 1; j < Ns; j++)//y
    {
        for (int i = 1; i < Ms; i++)//x
        {
            matrix[j][i] = get_aCoeff(Point(grid[j][j].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][j].x - 0.5 * h1, grid[j][i].y - 0.5 * h2), eps);
        }
    }
}

void get_bCoeffMatrix(const std::vector<std::vector<Point>>& grid,
        const double eps,
        std::vector<std::vector<double>>& matrix)
{
    int Ns = grid.size(), Ms = grid[0].size();
    double h1 = std::abs(grid[0][0].x - grid[0][1].x);
    double h2 = std::abs(grid[0][0].y - grid[1][0].y);
    for (int j = 1; j < Ns; j++)//y
    {
        for (int i = 1; i < Ms; i++)//x
        {
            matrix[j][i] = get_aCoeff(Point(grid[j][j].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][j].x + 0.5 * h1, grid[j][i].y + 0.5 * h2), eps);
        }
    }
}


void get_DiscMatrix(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& w,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2,
        std::vector<std::vector<double>>& matrix)
{
    // rk = Awk - B
    get_matrixA(aCoeff, bCoeff, w, h1, h2, matrix);
    int Ns = w.size(), Ms = w[0].size();
    for (int j = 1; j < Ns - 1; j++)//y
    {
        for (int i = 1; i < Ms - 1; i++)//x
        {
            matrix[j][i] -= B[j][i];
        }
    }
}




double get_IterParam(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2,
        const std::vector<std::vector<double>>& r_k,
        std::vector<std::vector<double>>& Ar)
{
    get_matrixA(aCoeff, bCoeff, r_k, h1, h2, Ar); // A примененное к rk
    return scal_prod(r_k, r_k, h1, h2) / scal_prod(Ar, r_k, h1, h2);
}