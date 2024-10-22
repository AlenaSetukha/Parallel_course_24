#ifndef _COMMON_TYPE_H_
#define _COMMON_TYPE_H_

#include <iostream>

struct Point
{
public:
    double x, y;
    Point();
    Point(const double x_in, const double y_in);
    Point(const Point& pnt_obj);
    ~Point(){};
};
#endif