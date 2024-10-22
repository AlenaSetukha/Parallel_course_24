#include "Point.hpp"

Point::Point()
{
    x = 0., y = 0.;
}

Point::Point(const double x_in, const double y_in)
{
    x = x_in, y = y_in;
}

Point::Point(const Point& pnt_obj)
{
    x = pnt_obj.x, y = pnt_obj.y;
}
