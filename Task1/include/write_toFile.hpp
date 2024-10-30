#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "Point.hpp"

void write_GridToFile(const std::string& result_dir, 
        const int Ms, const int Ns, 
        const std::vector<std::vector<Point>> grid);

void write_SolToFile(const std::string& result_dir, 
        const int Ms, const int Ns, 
        const std::vector<std::vector<double>> w_k);