#ifndef _COEFF_CALC_H_
#define _COEFF_CALC_H_

#include <iostream>
#include "Point.hpp"

double get_aCoeff(const Point& pnt1, const Point& pnt2, const double eps);

double get_bCoeff(const Point& pnt1, const Point& pnt2, const double eps);

double get_fCoeff(const double x1, const double x2, const double y1, const double y2);

#endif