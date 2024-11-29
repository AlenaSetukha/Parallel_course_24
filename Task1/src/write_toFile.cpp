#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "Point.hpp"

void write_GridToFile(const std::string& result_dir,
        const std::vector<std::vector<Point>> grid)
{
    std::ofstream fout_grid(result_dir + "grid.txt");
    if (!fout_grid.is_open()) {
        throw std::runtime_error("Open grid.txt error: main.cpp");
    }
    int Ns1 = grid.size(), Ms1 = grid[0].size();

    fout_grid << Ns1 << " " << Ms1 << std::endl;
    for (int j = 0; j < Ns1; j++) { //y
        for (int i = 0; i < Ms1; i++) { //x
            fout_grid << grid[j][i].x << " " << grid[j][i].y << std::endl;
        }
    }
    fout_grid.close();
}

void write_GridToFile(const std::string& result_dir, 
        const int Ms1, const int Ns1, const Point** grid)
{
    std::ofstream fout_grid(result_dir + "grid.txt");
    if (!fout_grid.is_open()) {
        throw std::runtime_error("Open grid.txt error: main.cpp");
    }
    fout_grid << Ns1 << " " << Ms1 << std::endl;
    for (int j = 0; j < Ns1; j++) { //y
        for (int i = 0; i < Ms1; i++) { //x
            fout_grid << grid[j][i].x << " " << grid[j][i].y << std::endl;
        }
    }
    fout_grid.close();
}



void write_SolToFile(const std::string& result_dir,
        const std::vector<std::vector<double>> w_k)
{
    std::ofstream fout_res(result_dir + "res.txt");
    if (!fout_res.is_open()) {
        throw std::runtime_error("Open res.txt error: main.cpp");
    }
    int Ns1 = w_k.size(), Ms1 = w_k[0].size();
    fout_res << Ns1 << " " << Ms1 << std::endl;
    for (int j = 0; j < Ns1; j++) { //y
        for (int i = 0; i < Ms1; i++) { //x
            fout_res << w_k[j][i] << std::endl;
        }
    }
    fout_res.close();
}








void write_MatrixToFile(const std::string& filename,
            const double** a_CoeffMatrix, const int Ns1, const int Ms1)
{
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Open aMatrix error: write_MatrixToFile.cpp");
    }
    fout << Ns1 << " " << Ms1 << std::endl;
    for (int j = 0; j < Ns1; j++) { //y
        for (int i = 0; i < Ms1; i++) { //x
            fout << a_CoeffMatrix[j][i] << std::endl;
        }
    }
    fout.close();
}

void write_MatrixToFile(const std::string& filename,
        const std::vector<std::vector<double>> a_CoeffMatrix)
{
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Open aMatrix error: write_MatrixToFile.cpp");
    }
    int Ns1 = a_CoeffMatrix.size(), Ms1 = a_CoeffMatrix[0].size();

    fout << Ns1 << " " << Ms1 << std::endl;
    for (int j = 0; j < Ns1; j++) { //y
        for (int i = 0; i < Ms1; i++) { //x
            fout << a_CoeffMatrix[j][i] << std::endl;
        }
    }
    fout.close();
}