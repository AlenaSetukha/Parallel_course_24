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
#include <mpi.h>
#include <time.h>

#include "Point.hpp"
#include "get_matrix.hpp"
#include "geom_calc.hpp"
#include "write_toFile.hpp"



/**
 * mpi + OpenMP решение задачи Дирихле для уравнения Пуассона в области D
 * методом фиктивных областей. Сделано через double**.
 * 
 * Количество доментов = количество MPI процессов
 * Кикл по k:
 *      1 шаг: Вычисление r[i][j] в каждом домене, вычисления независимы.
 * В конце шага нужно передать из процесса в процесс краевые значения,
 * которые лежат в общих частях доменов.
 *      2 шаг: вычисление тау в каждом домене: вычисляются числители и знаменатели.
 * В конце шага - сложить все числители и все знаменатели и разделить.
 *      3 шаг: вычисление W_k1 в каждом домене, независимо.
 *      4 шаг: получение нормы на каждом домене. В конце шага - найти максимум
 * среди всех норм на доменах. Если больше дельты, то идем снова на шаг 1.
 */




int main(int argc, char **argv)
{
    //===========Основные параметры задачи=======================
    const std::string param_fName = argc > 1 ? argv[1] : "../data/param_90_80.txt",\
                      result_dir = argc > 2 ? argv[2] : "../results_MPI/";
    int num_threads = argc > 3 ? std::stoi(argv[3]) : 4;
    auto start = std::chrono::high_resolution_clock::now();
    int Ns, Ms, k_max;
    double A1, B1, A2, B2, delta;

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    //=====================Считывание параметров==================
    std::ifstream fin(param_fName);
    if (!fin.is_open()) {
        throw std::runtime_error("Open param_fName error: main.cpp");
    }
    fin >> Ns >> Ms >> A1 >> B1 >> A2 >> B2 >> delta >> k_max;
    fin.close();

    if (Ns % 2 != 0 || Ms % 2 != 0) {
        throw std::runtime_error("Число ячеек MxN должно быть четным!");
    }


    //==========================Сетка=============================
    const double h1 = (B1 - A1) / Ms;
    const double h2 = (B2 - A2) / Ns;
    const double eps = std::max(h1, h2) * std::max(h1, h2);
    Point** grid = new Point*[Ns + 1];
    for (int j = 0; j < Ns + 1; j++)
    {
        grid[j] = new Point[Ms + 1];
        for (int i = 0; i < Ms + 1; i++) { //x
            grid[j][i] = Point(A1 + i * h1, B2 - j * h2);
        }
    }


    //=================Набор константных матриц===================
    //==========Обязательно чтобы выделение было подряд===========
    double** B = new double*[Ns + 1];
    for (int j = 0; j < Ns + 1; j++) {
        B[j] = new double[Ms + 1];
        std::fill(B[j], B[j] + Ms + 1, 0.0);
    }

    double** a_CoeffMatrix = new double*[Ns + 1];
    for (int j = 0; j < Ns + 1; j++) {
        a_CoeffMatrix[j] = new double[Ms + 1];
        std::fill(a_CoeffMatrix[j], a_CoeffMatrix[j] + Ms + 1, 0.0);
    }

    double** b_CoeffMatrix =  new double*[Ns + 1];
    for (int j = 0; j < Ns + 1; j++) {
        b_CoeffMatrix[j] = new double[Ms + 1];
        std::fill(b_CoeffMatrix[j], b_CoeffMatrix[j] + Ms + 1, 0.0);
    }

    get_constMatrixOMP((const Point**)grid, eps, Ns + 1, Ms + 1, B, a_CoeffMatrix, b_CoeffMatrix);




    //======================Разделение на процессы==================
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);     // Получаем количество MPI процессов
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);     // Номер текущего процесса

    if (!(size == 2 || size == 4))
    {
        // Заверение всех MPI процессов
        MPI_Finalize();   
        throw std::runtime_error("Число MPI процессов не 2 и не 4!");
    }


    //===============Определение текущего поддомена=================
    // Число узлов по x и y в домене, нумерация с нуля(так же, как и Ms, Ns)
    int M_ofDomen = std::ceil((Ms + 1) / size) + 1;
    int N_ofDomen = Ns; 

    int startX_Indx, startY_Indx; //  (0,0) либо (3, 0)
    if (size == 2) {
        startX_Indx = (rank % 2) * ((Ms / 2) - 1);
        startY_Indx = 0;
    }


    std::cout << "Номер домена: " << rank << " , N и M: " << N_ofDomen << " " << M_ofDomen << std::endl;
    std::cout << " Индексы начала домена по y и по x:  " << startY_Indx << " " << startX_Indx << std::endl;
    

    //====================Итерационный процесс=======================
    double** w_k_loc = new double*[N_ofDomen + 1];
    double** w_k1_loc = new double*[N_ofDomen + 1];
    double** r_k_loc = new double*[N_ofDomen + 1];
    double** Ar_loc = new double*[N_ofDomen + 1];
    for (int j = 0; j < N_ofDomen + 1; j++) {
        w_k_loc[j] = new double[M_ofDomen + 1];  // стартовое приближение
        w_k1_loc[j] = new double[M_ofDomen + 1]; // приближение на (k+1) шаге
        r_k_loc[j] = new double[M_ofDomen + 1];  // невязка Awk - B
        Ar_loc[j] = new double[M_ofDomen + 1];   // Ar - заполняется на каждом шаге
        std::fill(w_k_loc[j], w_k_loc[j] + M_ofDomen + 1, 0.0);
        std::fill(w_k1_loc[j], w_k1_loc[j] + M_ofDomen + 1, 0.0);
        std::fill(r_k_loc[j], r_k_loc[j] + M_ofDomen + 1, 0.0);
        std::fill(Ar_loc[j], Ar_loc[j] + M_ofDomen + 1, 0.0);
    }


    double* buf_r = new double[N_ofDomen + 1];
    double* buf_w1 = new double[N_ofDomen + 1];

    double global_norm;
    const int _TAG_R01 = 15;
    const int _TAG_R10 = 16;
    const int _TAG_W01 = 25;
    const int _TAG_W10 = 26;

    for (int k = 1; k < k_max; k++)
    {
        /**-------------------------------------------------------------
         *  Заполнение локальной невязки на домене rk = Awk - B.
         *  Заполнение происходит во всех внутренних узлах области,
         *  но с использованием w_k во всех узлах области(крайних тоже).
         *-------------------------------------------------------------*/
        get_rkOMP((const double**)(a_CoeffMatrix + startY_Indx),
                (const double**)(b_CoeffMatrix + startY_Indx),
                (const double**)(B + startY_Indx), startX_Indx, (const double**)w_k_loc,
                N_ofDomen, M_ofDomen, h1, h2, r_k_loc);
        

        /**-------------------------------------------------------------
         * ОБМЕН МЕЖДУ ПРОЦЕССАМИ r[i][j] для определния rk_loc[i][j]
         * на правой или левой границе(на всех остальных границах = ноль)
         *-------------------------------------------------------------*/
        MPI_Status status;
        if (!rank) {
            for (int j = 0; j < N_ofDomen + 1; j++) {
                buf_r[j] = r_k_loc[j][M_ofDomen - 2]; //записали левый край для второго домена
            }
            MPI_Send(&buf_r[0], N_ofDomen + 1, MPI_DOUBLE, 1, _TAG_R01, MPI_COMM_WORLD);
            MPI_Recv(&buf_r[0], N_ofDomen + 1, MPI_DOUBLE, 1, _TAG_R10, MPI_COMM_WORLD, &status);
            for (int j = 0; j < N_ofDomen + 1; j++) {
                r_k_loc[j][M_ofDomen] = buf_r[j]; //переложили правый край домена
            }
        } else if (rank == 1) {
            for (int j = 0; j < N_ofDomen + 1; j++) {
                buf_r[j] = r_k_loc[j][2];           //записали правый край для первого домена
            }
            MPI_Send(&buf_r[0], N_ofDomen + 1, MPI_DOUBLE, 0, _TAG_R10, MPI_COMM_WORLD);    
            MPI_Recv(&buf_r[0], N_ofDomen + 1, MPI_DOUBLE, 0, _TAG_R01, MPI_COMM_WORLD, &status);
            for (int j = 0; j < N_ofDomen + 1; j++) {
                r_k_loc[j][0] = buf_r[j]; //переложили левый край домена
            }
        }
        /**-------------------------------------------------------------
         * Вычисление числителя и знаменателя итерационного параметра в
         * текущем домене:
         *      tau_k1 = (rk, rk) / (Ark, rk) = num / denom
         * Скалярное произведение считается во всех внутренних узлах,
         * оператор Ar - во внутренних, но использованием всех. 
         *-------------------------------------------------------------*/
        double tau_loc_num = 0., tau_loc_denom = 0.;
        
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < N_ofDomen; j++) {
            for (int i = 1; i < M_ofDomen; i++) {
                double diff1 = a_CoeffMatrix[j + startY_Indx][startX_Indx + i + 1] *
                    (r_k_loc[j][i + 1] - r_k_loc[j][i]) -
                    a_CoeffMatrix[j + startY_Indx][startX_Indx + i] *
                    (r_k_loc[j][i] - r_k_loc[j][i - 1]);// aCoeff[j][i + 1] * (wk[j][i + 1] - wk[j][i]) - aCoeff[j][i] * (wk[j][i] - wk[j][i - 1]);
                double diff2 = b_CoeffMatrix[j + startY_Indx +  1][startX_Indx + i]*
                    (r_k_loc[j + 1][i] - r_k_loc[j][i]) -
                    b_CoeffMatrix[j + startY_Indx][startX_Indx + i] *
                    (r_k_loc[j][i] - r_k_loc[j - 1][i]);
                Ar_loc[j][i] = -diff1 / (h1 * h1) - diff2 / (h2 * h2);
            }
        }


        if (!rank) {
            for (int j = 1; j < N_ofDomen; j++) {
                for (int i = 1; i < M_ofDomen - 1; i++) {
                    tau_loc_denom += Ar_loc[j][i] * r_k_loc[j][i];
                    tau_loc_num += r_k_loc[j][i] * r_k_loc[j][i];
                }
            }
        } else {
            for (int j = 1; j < N_ofDomen; j++) {
                for (int i = 1; i < M_ofDomen; i++) {
                    tau_loc_denom += Ar_loc[j][i] * r_k_loc[j][i];
                    tau_loc_num += r_k_loc[j][i] * r_k_loc[j][i];
                }
            }
        }
        tau_loc_num *= h1 * h2;
        tau_loc_denom *= h1 * h2;

        /**-------------------------------------------------------------
         * СЛИЯНИЕ ПРОЦЕССОВ для получения общего итерационного параметра.
         *-------------------------------------------------------------*/
        double global_num, global_denom;
        MPI_Allreduce(&tau_loc_num, &global_num, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&tau_loc_denom, &global_denom, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double global_tau = global_num / global_denom;

        /**-------------------------------------------------------------
         * Сдвиг решения на домене во внутренних узлах:
         *      wk1 = wk - tau_k1 * rk
         *-------------------------------------------------------------*/
        for (int j = 1; j < N_ofDomen; j++) { //y
            for (int i = 1; i < M_ofDomen; i++) { //x
                w_k1_loc[j][i] = w_k_loc[j][i] - global_tau * r_k_loc[j][i];
            }
        }



        /**-------------------------------------------------------------
         * ОБМЕН МЕЖДУ ПРОЦЕССАМИ для определния wk1_loc[i][j] на правой
         * или левой границе(на всех остальных границах = ноль)
         *-------------------------------------------------------------*/
        if (!rank) {
            for (int j = 0; j < N_ofDomen + 1; j++) {
                buf_w1[j] = w_k1_loc[j][M_ofDomen - 2]; //записали левый край для второго домена
            }
            MPI_Send(&buf_w1[0], N_ofDomen + 1, MPI_DOUBLE, 1, _TAG_W01, MPI_COMM_WORLD);
            MPI_Recv(&buf_w1[0], N_ofDomen + 1, MPI_DOUBLE, 1, _TAG_W10, MPI_COMM_WORLD, &status);
            for (int j = 0; j < N_ofDomen + 1; j++) {
                w_k1_loc[j][M_ofDomen] = buf_w1[j]; //переложили правый край домена
            }
        } else if (rank == 1) {
            for (int j = 0; j < N_ofDomen + 1; j++) {
                buf_w1[j] = w_k1_loc[j][2];           //записали правый край для первого домена
            }
            MPI_Send(&buf_w1[0], N_ofDomen + 1, MPI_DOUBLE, 0, _TAG_W10, MPI_COMM_WORLD);      
            MPI_Recv(&buf_w1[0], N_ofDomen + 1, MPI_DOUBLE, 0, _TAG_W01, MPI_COMM_WORLD, &status);
            for (int j = 0; j < N_ofDomen + 1; j++) {
                w_k1_loc[j][0] = buf_w1[j]; //переложили левый край домена
            }     
        }
        /**-------------------------------------------------------------
         * Вычисление локальной нормы и СЛИЯНИЕ ПРОЦЕССОВ для получения
         * общей нормы решения. 
         *-------------------------------------------------------------*/
        double norm_loc = get_normC((const double**)w_k1_loc, (const double**)w_k_loc, N_ofDomen, M_ofDomen);
        MPI_Allreduce(&norm_loc, &global_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);



        if (global_norm < delta)
        {
            std::cout << "Достигнута точность! Количество шагов: " << k << std::endl;
            break;
        } else {
            // Перекладываем все значения
            for (int j = 0; j < N_ofDomen + 1; j++) { //y
                for (int i = 0; i < M_ofDomen + 1; i++) { //x
                    w_k_loc[j][i] = w_k1_loc[j][i];
                }
            }
        }
    }

    //======Запись локального результата====
    std::ofstream fout_res(result_dir + "res" + std::to_string(rank) + ".txt");
    std::ofstream fout_gridLoc(result_dir + "grid" + std::to_string(rank) + ".txt");
    if (!fout_res.is_open() && !!fout_gridLoc.is_open()) {
        throw std::runtime_error("Open res.txt error: main.cpp");
    }


    if (!rank) {
        fout_res << N_ofDomen + 1 << " " << M_ofDomen << std::endl;
        fout_gridLoc << N_ofDomen + 1 << " " << M_ofDomen << std::endl;
        for (int j = 0; j < N_ofDomen + 1; j++) { //y
            for (int i = 0; i < M_ofDomen; i++) { //x
                fout_res << w_k1_loc[j][i] << std::endl;
                fout_gridLoc << grid[j + startY_Indx][i + startX_Indx].x << " " << grid[j + startY_Indx][i + startX_Indx].y << std::endl;
            }
        }
    } else if (rank == 1) {
        fout_res << N_ofDomen + 1 << " " << M_ofDomen << std::endl;
        fout_gridLoc << N_ofDomen + 1 << " " << M_ofDomen << std::endl;
        for (int j = 0; j < N_ofDomen + 1; j++) { //y
            for (int i = 1; i < M_ofDomen + 1; i++) { //x
                fout_res << w_k1_loc[j][i] << std::endl;
                fout_gridLoc << grid[j + startY_Indx][i + startX_Indx].x << " " << grid[j + startY_Indx][i + startX_Indx].y << std::endl;
            }
        }
    }
    fout_res.close();
    fout_gridLoc.close();


    // Очистка памяти всех процессов
    for (int j = 0; j < N_ofDomen + 1; j++)
    {
        delete[] w_k_loc[j];  // стартовое приближение
        delete[] w_k1_loc[j]; // приближение на (k+1) шаге
        delete[] r_k_loc[j];  // невязка Awk - B
        delete[] Ar_loc[j];   // Ar - заполняется на каждом шаге
    }
    delete[] w_k1_loc;
    delete[] w_k_loc;
    delete[] r_k_loc;
    delete[] Ar_loc;

    delete[] buf_r;
    delete[] buf_w1;

    //========Запись ответа========
    if (!rank) {
        std::cout << "Предельное число шагов: " << k_max << std::endl;
        std::cout << "Достигнутая норма: " << global_norm << std::endl;
        //write_GridToFile(result_dir, Ms + 1, Ns + 1, (const Point**)grid);

        //========Расчет времени=======
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    }



    // Очистка памяти
    for (int j = 0; j < Ns + 1; j++)
    {
        delete[] grid[j];
        delete[] B[j];
        delete[] a_CoeffMatrix[j];
        delete[] b_CoeffMatrix[j];
    }
    delete[] grid;
    delete[] B;
    delete[] a_CoeffMatrix;
    delete[] b_CoeffMatrix;

    // Завершение всех MPI процессов
    MPI_Finalize();
    return 0;
}