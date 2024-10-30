#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "Point.hpp"

void write_GridToFile(const std::string& result_dir, 
        const int Ms, const int Ns, 
        const std::vector<std::vector<Point>> grid)
{
    std::ofstream fout_grid(result_dir + "grid.txt");
    if (!fout_grid.is_open()) {
        throw std::runtime_error("Open grid.txt error: main.cpp");
    }
    fout_grid << Ns - 1 << " " << Ms - 1 << std::endl;
    for (int j = 1; j < Ns; j++) { //y
        for (int i = 1; i < Ms; i++) { //x
            fout_grid << grid[j][i].x << " " << grid[j][i].y << std::endl;
        }
    }
    fout_grid.close();
}

void write_SolToFile(const std::string& result_dir, 
        const int Ms, const int Ns, 
        const std::vector<std::vector<double>> w_k)
{
    std::ofstream fout_res(result_dir + "res.txt");
    if (!fout_res.is_open()) {
        throw std::runtime_error("Open res.txt error: main.cpp");
    }
    fout_res << Ns - 1 << " " << Ms - 1 << std::endl;
    for (int j = 1; j < Ns; j++) { //y
        for (int i = 1; i < Ms; i++) { //x
            fout_res << w_k[j][i] << std::endl;
        }
    }
    fout_res.close();
}