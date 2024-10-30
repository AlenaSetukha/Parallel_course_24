#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>
#include <string>

#include "Point.hpp"
#include "get_matrix.hpp"
#include "geom_calc.hpp"
#include "write_toFile.hpp"

/**
 * Последовательное решение задачи Дирихле для уравнения Пуассона в области D
 * методом фиктивных областей.
 * Задача решается в фиктивном прямоугольнике
 *      [A1 ; B1] x [A2 ; B2] = [-1.5; 1.5] x [-1; 1].
 * 
 * Область D задается уравнением эллипса:
 *      x^2 + 4y^2 = 1
 * 
 * Разностная схема - метод конечных разностей.
 * Начало сетки - правый верхний угол.
 */

int main(int argc, char **argv)
{
    const std::string param_fName = argc > 1 ? argv[1] : "../data/param_40.txt",\
                      result_dir = argc > 2 ? argv[2] : "../results/";
    auto start = std::chrono::high_resolution_clock::now();

    int Ns, Ms, k_max;
    double A1, B1, A2, B2, delta;

    //====Считывание параметров====
    std::ifstream fin(param_fName);
    if (!fin.is_open()) {
        throw std::runtime_error("Open param_fName error: main.cpp");
    }
    fin >> Ns >> Ms >> A1 >> B1 >> A2 >> B2 >> delta >> k_max;
    fin.close();
    
    //===========Сетка=============
    const double h1 = std::abs(B1 - A1) / Ms;
    const double h2 = std::abs(B2 - A2) / Ns;
    const double eps = std::max(h1, h2) * std::max(h1, h2);
    std::vector<std::vector<Point>> grid(Ns + 1, std::vector<Point>(Ms + 1));
    for (int j = 0; j < Ns + 1; j++) { //y
        for (int i = 0; i < Ms + 1; i++) { //x
            grid[j][i] = Point(A1 + i * h1, B2 - j * h2);
        }
    }

    //========Запись сетки=========
    write_GridToFile(result_dir, Ms, Ns, grid);



    //===Набор константных матриц==
    std::vector<std::vector<double>> B(Ns + 1, std::vector<double>(Ms + 1, 0.));          // матрица B = F[i][j]
    get_matrixB(grid, B);

    std::vector<std::vector<double>> a_CoeffMatrix(Ns + 1, std::vector<double>(Ms + 1, 0.));  // матрица a[i][j]
    std::vector<std::vector<double>> b_CoeffMatrix(Ns + 1, std::vector<double>(Ms + 1, 0.));  // матрица b[i][j]
    get_aCoeffMatrix(grid, eps, a_CoeffMatrix);
    get_bCoeffMatrix(grid, eps, b_CoeffMatrix);

    //=====Итерационный процесс====
    std::vector<std::vector<double>> w_k(Ns + 1, std::vector<double>(Ms + 1, 0.));      // стартовое приближение
    std::vector<std::vector<double>> w_k1(Ns + 1, std::vector<double>(Ms + 1, 0.));     // приближение на (k+1) шаге
    std::vector<std::vector<double>> r_k(Ns + 1, std::vector<double>(Ms + 1, 0.));      // невязка Awk - B
    std::vector<std::vector<double>> Ar(Ns + 1, std::vector<double>(Ms + 1, 0.));       // Ar - заполняется на каждом шаге заново

    double tau_k1, norm;
    for (int k = 1; k < k_max; k++)
    {
        get_DiscMatrix(a_CoeffMatrix, b_CoeffMatrix, w_k, B,
                h1, h2, r_k);
        tau_k1 = get_IterParam(a_CoeffMatrix, b_CoeffMatrix, B,
                h1, h2, r_k, Ar);
        for (int j = 1; j < Ns - 1; j++) { //y
            for (int i = 1; i < Ms - 1; i++) { //x
                w_k1[j][i] = w_k[j][i] - tau_k1 * r_k[j][i];
            }
        }
        norm = get_normC(w_k1, w_k); //сделать поиск максимум по всем узлам, не только внутренним
if (k == 10)
{
    std::cout << "На 10-м шаге: " << std::endl;
    std::cout << "      w_10[2][5] = " << w_k1[2][5] << std::endl;
    std::cout << "      w_10[45][40] = " << w_k1[45][40] << std::endl;
    std::cout << "      w_10[50][40] = " << w_k1[50][40] << std::endl;
    std::cout << "      norm = " << norm << std::endl;
}
        if (norm < delta)
        {
            std::cout << "Достигнута точность! Количество шагов: " << k << std::endl;
            break;
        } else {
            for (int j = 1; j < Ns - 1; j++) { //y
                for (int i = 1; i < Ms - 1; i++) { //x
                    w_k[j][i] = w_k1[j][i];
                }
            }
        }
    }
    //========Запись ответа========
    write_SolToFile(result_dir, Ms, Ns, w_k);
    std::cout << "Предельное число шагов: " << k_max << std::endl;
    std::cout << "Достигнутая норма: " << norm << std::endl;


    //========Расчет времени=======
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}
