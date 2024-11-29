#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <omp.h> 

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
            matrix[j][i] = get_aCoeff(Point(grid[j][i].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][i].x - 0.5 * h1, grid[j][i].y - 0.5 * h2), eps);
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
            matrix[j][i] = get_bCoeff(Point(grid[j][i].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][i].x + 0.5 * h1, grid[j][i].y + 0.5 * h2), eps);
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










//===============Параллельный набор матриц=========================
// Применение оператора A к вектору wk во всех внутренних узлах
void get_matrixA_OMP(const double** aCoeff, const double** bCoeff,
        const int start_XIndx, const double** wk, const double h1, const double h2,
        const int N_ofDomen, const int M_ofDomen, double** r_k)
{
    // rk = Awk
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < N_ofDomen; j++)//y
    {
        for (int i = 1; i < M_ofDomen; i++)//x
        {
            double diff1 = aCoeff[j][start_XIndx + i + 1] *
                        (wk[j][i + 1] - wk[j][i]) -
                        aCoeff[j][start_XIndx + i] *
                        (wk[j][i] - wk[j][i - 1]);// aCoeff[j][i + 1] * (wk[j][i + 1] - wk[j][i]) - aCoeff[j][i] * (wk[j][i] - wk[j][i - 1]);
            double diff2 = bCoeff[j + 1][start_XIndx + i]*
                        (wk[j + 1][i] - wk[j][i]) -
                        bCoeff[j][start_XIndx + i] *
                        (wk[j][i] - wk[j - 1][i]);
            r_k[j][i] = -diff1 / sqr(h1) - diff2 / sqr(h2);
        }
    }
}




// Для векторов
void get_constMatrixOMP(const std::vector<std::vector<Point>>& grid,
        const double eps, std::vector<std::vector<double>>& B,
        std::vector<std::vector<double>>& a_CoeffMatrix,
        std::vector<std::vector<double>>& b_CoeffMatrix)
{
    int Ns1 = grid.size(), Ms1 = grid[0].size();
    double h1 = std::abs(grid[0][0].x - grid[0][1].x);
    double h2 = std::abs(grid[0][0].y - grid[1][0].y);
    #pragma omp parallel for
    for (int j = 1; j < Ns1; j++)//y
    {
        for (int i = 1; i < Ms1; i++)//x
        {
            a_CoeffMatrix[j][i] = get_aCoeff(Point(grid[j][i].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][i].x - 0.5 * h1, grid[j][i].y - 0.5 * h2), eps);
            b_CoeffMatrix[j][i] = get_bCoeff(Point(grid[j][i].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][i].x + 0.5 * h1, grid[j][i].y + 0.5 * h2), eps);
            if (j != Ns1 - 1 && i != Ms1 - 1)
            {
                B[j][i] = get_fCoeff(grid[j][i].x - 0.5 * h1, grid[j][i].x + 0.5 * h1,
                                      grid[j][i].y + 0.5 * h2, grid[j][i].y - 0.5 * h2);
            }
        }
    }
}


// Для указателей double**
void get_constMatrixOMP(const Point** grid, const double eps,
        const int Ns1, const int Ms1, double** B,
        double** a_CoeffMatrix, double** b_CoeffMatrix)
{
    double h1 = std::abs(grid[0][0].x - grid[0][1].x);
    double h2 = std::abs(grid[0][0].y - grid[1][0].y);
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < Ns1; j++)//y
    {
        for (int i = 1; i < Ms1; i++)//x
        {
            a_CoeffMatrix[j][i] = get_aCoeff(Point(grid[j][i].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][i].x - 0.5 * h1, grid[j][i].y - 0.5 * h2), eps);
            b_CoeffMatrix[j][i] = get_bCoeff(Point(grid[j][i].x - 0.5 * h1, grid[j][i].y + 0.5 * h2),
                                      Point(grid[j][i].x + 0.5 * h1, grid[j][i].y + 0.5 * h2), eps);
            if (j != Ns1 - 1 && i != Ms1 - 1)
            {
                B[j][i] = get_fCoeff(grid[j][i].x - 0.5 * h1, grid[j][i].x + 0.5 * h1,
                                      grid[j][i].y + 0.5 * h2, grid[j][i].y - 0.5 * h2);
            }
        }
    }
}





// Для векторов
void get_rkOMP(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& wk,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2,
        std::vector<std::vector<double>>& r_k)
{
    // rk = Awk - B
    int Ns1 = wk.size(), Ms1 = wk[0].size();
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < Ns1 - 1; j++)//y
    {
        for (int i = 1; i < Ms1 - 1; i++)//x
        {
            double diff1 = aCoeff[j][i + 1] * (wk[j][i + 1] - wk[j][i]) -
                        aCoeff[j][i] * (wk[j][i] - wk[j][i - 1]);
            double diff2 = bCoeff[j + 1][i] * (wk[j + 1][i] - wk[j][i]) -
                        bCoeff[j][i] * (wk[j][i] - wk[j - 1][i]);
            r_k[j][i] = -diff1 / sqr(h1) - diff2 / sqr(h2) - B[j][i];
        }
    }
}


// Для  указателей. Матрица в виде ленты.
// Ширина матрицы: M_ofDomen + 1, длина матрицы: N_ofDomen + 1
// Ns, Ms - глубина и шрина общих матриц
void get_rkOMP(const double** aCoeff, const double** bCoeff,
        const double** B, const int start_Xindx, const double** wk,
        const int N_ofDomen, const int M_ofDomen,
        const double h1, const double h2, double** r_k)
{
    // rk = Awk - B
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < N_ofDomen; j++)//y
    {
        for (int i = 1; i < M_ofDomen; i++)//x
        {
            double diff1 = aCoeff[j][start_Xindx + i + 1] *
                        (wk[j][i + 1] - wk[j][i]) -
                        aCoeff[j][start_Xindx + i] *
                        (wk[j][i] - wk[j][i - 1]);  // aCoeff[j][i + 1] * (wk[j][i + 1] - wk[j][i]) - aCoeff[j][i] * (wk[j][i] - wk[j][i - 1]);
            double diff2 = bCoeff[j + 1][start_Xindx + i] *
                        (wk[j + 1][i] - wk[j][i]) -
                        bCoeff[j][start_Xindx + i] *
                        (wk[j][i] - wk[j - 1][i]);

            r_k[j][i] = -diff1 / sqr(h1) - diff2 / sqr(h2) - B[j][start_Xindx + i];
        }
    }
}









double get_IterParam_OMP(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const double h1, const double h2,
        const std::vector<std::vector<double>>& r_k,
        std::vector<std::vector<double>>& Ar)
{
    // A примененное к rk
    int Ns1 = r_k.size(), Ms1 = r_k[0].size();
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < Ns1 - 1; j++)//y
    {
        for (int i = 1; i < Ms1 - 1; i++)//x
        {
            double diff1 = aCoeff[j][i + 1] * (r_k[j][i + 1] - r_k[j][i]) -
                        aCoeff[j][i] * (r_k[j][i] - r_k[j][i - 1]);
            double diff2 = bCoeff[j + 1][i] * (r_k[j + 1][i] - r_k[j][i]) -
                        bCoeff[j][i] * (r_k[j][i] - r_k[j - 1][i]);
            Ar[j][i] = -diff1 / sqr(h1) - diff2 / sqr(h2);
        }
    }
    return scal_prod(r_k, r_k, h1, h2) /
                 scal_prod(Ar, r_k, h1, h2);
}






void sol_StepOMP(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& wk,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2,
        std::vector<std::vector<double>>& r_k,
        std::vector<std::vector<double>>& Ar,
        std::vector<std::vector<double>>& wk1)
{
    // rk = Awk - B
    get_rkOMP(aCoeff, bCoeff, wk, B, h1, h2, r_k);
    // tau_k1 = (rk, rk) / (Ark, rk)
    double tau_k1 = get_IterParam_OMP(aCoeff, bCoeff, h1, h2, r_k, Ar);

std::cout << tau_k1 << std::endl;

    // wk1 = wk - tau_k1 * rk
    int Ns1 = wk.size(), Ms1 = wk[0].size();
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < Ns1; j++) { //y
        for (int i = 0; i < Ms1; i++) { //x
            wk1[j][i] = wk[j][i] - tau_k1 * r_k[j][i];
        }
    }
}
