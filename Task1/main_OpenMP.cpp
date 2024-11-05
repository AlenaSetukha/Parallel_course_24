#include <omp.h> 
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>
#include <string>
#include <stdexcept>

#include "Point.hpp"
#include "get_matrix.hpp"
#include "geom_calc.hpp"
#include "write_toFile.hpp"

/**
 * OpenMP решение задачи Дирихле для уравнения Пуассона в области D
 * методом фиктивных областей.
 * Распараллелены по потокам:
 *      - набор сетки
 *      - набор статических матриц: B(=F), a, b
 *      - набор матрицы A на каждом шаге
 *      - скалярное произведение
 *      - поиск C-нормы
 */

int main(int argc, char **argv)
{
    const std::string param_fName = argc > 1 ? argv[1] : "../data/param.txt",\
                      result_dir = argc > 2 ? argv[2] : "../results/";
    int num_threads = argc > 3 ? std::stoi(argv[3]) : 4;
    auto start = std::chrono::high_resolution_clock::now();
    int Ns, Ms, k_max;
    double A1, B1, A2, B2, delta;

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    //====Считывание параметров====
    std::ifstream fin(param_fName);
    if (!fin.is_open()) {
        throw std::runtime_error("Open param_fName error: main.cpp");
    }
    fin >> Ns >> Ms >> A1 >> B1 >> A2 >> B2 >> delta >> k_max;
    fin.close();


    //===========Сетка=============
    const double h1 = (B1 - A1) / Ms;
    const double h2 = (B2 - A2) / Ns;
    const double eps = std::max(h1, h2) * std::max(h1, h2);
    std::vector<std::vector<Point>> grid(Ns + 1, std::vector<Point>(Ms + 1));

    #pragma omp parallel for
    for (int j = 0; j < Ns + 1; j++) { //y
        for (int i = 0; i < Ms + 1; i++) { //x
            grid[j][i] = Point(A1 + i * h1, B2 - j * h2);
        }
    }
    write_GridToFile(result_dir, Ms, Ns, grid);




    //===Набор константных матриц==
    std::vector<std::vector<double>> B(Ns + 1, std::vector<double>(Ms + 1, 0.));          // матрица B = F[i][j]
    std::vector<std::vector<double>> a_CoeffMatrix(Ns + 1, std::vector<double>(Ms + 1, 0.));  // матрица a[i][j]
    std::vector<std::vector<double>> b_CoeffMatrix(Ns + 1, std::vector<double>(Ms + 1, 0.));  // матрица b[i][j]
    get_constMatrixOMP(grid, eps, B, a_CoeffMatrix, b_CoeffMatrix);

    //=====Итерационный процесс====
    std::vector<std::vector<double>> w_k(Ns + 1, std::vector<double>(Ms + 1, 0.));      // стартовое приближение
    std::vector<std::vector<double>> w_k1(Ns + 1, std::vector<double>(Ms + 1, 0.));     // приближение на (k+1) шаге
    std::vector<std::vector<double>> r_k(Ns + 1, std::vector<double>(Ms + 1, 0.));      // невязка Awk - B
    std::vector<std::vector<double>> Ar(Ns + 1, std::vector<double>(Ms + 1, 0.));       // Ar - заполняется на каждом шаге заново

    double norm;
    for (int k = 1; k < k_max; k++)
    {
        sol_StepOMP(a_CoeffMatrix, b_CoeffMatrix, w_k, B,
                h1, h2, r_k, Ar, w_k1);

        norm = get_normC_OMP(w_k1, w_k);

        if (norm < delta)
        {
            std::cout << "Достигнута точность! Количество шагов: " << k << std::endl;
            break;
        } else {
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < Ns - 1; j++) { //y
                for (int i = 1; i < Ms - 1; i++) { //x
                    w_k[j][i] = w_k1[j][i];
                }
            }
        }
    }
    std::cout << "Предельное число шагов: " << k_max << std::endl;
    std::cout << "Достигнутая норма: " << norm << std::endl;
    //========Запись ответа========
    write_SolToFile(result_dir, Ms, Ns, w_k);

    //========Расчет времени=======
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}