#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "Point.hpp"

void write_GridToFile(const std::string& result_dir,
        const std::vector<std::vector<Point>> grid);

void write_GridToFile(const std::string& result_dir, 
        const int Ms1, const int Ns1, const Point** grid);

void write_SolToFile(const std::string& result_dir,
        const std::vector<std::vector<double>> w_k);


void write_MatrixToFile(const std::string& filename,
        const double** a_CoeffMatrix, const int Ns1, const int Ms1);

void write_MatrixToFile(const std::string& filename,
        const std::vector<std::vector<double>> a_CoeffMatrix);