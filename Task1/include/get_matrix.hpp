#ifndef _GET_MATRIX_H_
#define _GET_MATRIX_H_

#include <complex>
#include <vector>
#include "Point.hpp"

void get_matrixA(const std::vector<std::vector<double>>& a,
        const std::vector<std::vector<double>>& b,
        const std::vector<std::vector<double>>& w,
        const double h1, const double h2,
        std::vector<std::vector<double>>& matrix);





//======================Константные матрицы==========================
//====================последовательный набор=========================
void get_matrixB(const std::vector<std::vector<Point>>& grid,
        std::vector<std::vector<double>>& matrix);


void get_aCoeffMatrix(const std::vector<std::vector<Point>>& grid,
        const double eps,
        std::vector<std::vector<double>>& matrix);


void get_bCoeffMatrix(const std::vector<std::vector<Point>>& grid,
        const double eps,
        std::vector<std::vector<double>>& matrix);





void get_DiscMatrix(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& w,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2,
        std::vector<std::vector<double>>& matrix);


double get_IterParam(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2,
        const std::vector<std::vector<double>>& r_k,
        std::vector<std::vector<double>>& Ar);
        














//======================Константные матрицы==========================
//=======================параллельный набор==========================
void get_constMatrixOMP(const std::vector<std::vector<Point>>& grid,
        const double eps, const int num_threads, 
        std::vector<std::vector<double>>& B,
        std::vector<std::vector<double>>& a_CoeffMatrix,
        std::vector<std::vector<double>>& b_CoeffMatrix);



void get_rkOMP(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& wk,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2, const int num_threads,
        std::vector<std::vector<double>>& r_k);


double get_IterParam_OMP(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const double h1, const double h2,
        const std::vector<std::vector<double>>& r_k,
        const int num_threads, std::vector<std::vector<double>>& Ar);

void sol_StepOMP(const std::vector<std::vector<double>>& aCoeff,
        const std::vector<std::vector<double>>& bCoeff,
        const std::vector<std::vector<double>>& w,
        const std::vector<std::vector<double>>& B,
        const double h1, const double h2, const int num_threads,
        std::vector<std::vector<double>>& r_k,
        std::vector<std::vector<double>>& Ar,
        std::vector<std::vector<double>>& wk1);

#endif