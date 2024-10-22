#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "Point.hpp"
#include "geom_calc.hpp"


static inline double sqr(const double x)
{
    return x * x;
}

static bool inEllips(const double x, const double y)
{
    return sqr(x) + 4.0 * sqr(y) <= 1.;
}

/**
 * Функция, определяющая коэффициент a[i][j] по положению
 * отрезка относительно кривой(он гарантированно состоит из
 * внутренних узлов сетки).
 * Уравнение эллипса:
 *      x^2 + 4y^2 = 1
 * Внешняя фиктивная область:
 *      (x,y) = [-1.5; 1.5] x [-1; 1]
 * Возвращает:
 *      1 - если отрезок полностью в области D
 *      1 / eps - если отрезок полностью в области D^
 *      l/h2 + (1 - l/h2) / eps - если отрезок пересекает Gamma
 */
double get_aCoeff(const Point& pnt1, const Point& pnt2, const double eps)
{
    bool p1 = inEllips(pnt1.x, pnt1.y), p2 = inEllips(pnt2.x, pnt2.y);

    // Обе точки внутри эллипса
    if (p1  && p2) {
        return 1.;
    } 
    
    // Проверка наличия точек пересечения с эллипсом
    double deg = (1.0 - sqr(pnt1.x)) / 4.0;
    if (deg < 1e-9) { // Обе точки снаружи, без пересечений
        return 1.0 / eps;
    }
    
    double y3p = sqrt(deg), y3m = -sqrt(deg);
    double h2 = std::abs(pnt2.y - pnt1.y);
    double min_y = std::min(pnt2.y, pnt1.y), max_y = std::max(pnt2.y, pnt1.y);
    
    
    // Обе точки снаружи
    if (!p1 && !p2) {
        if (y3p <= min_y || y3m >= max_y) {
            // Без пересечений
            return 1.0 / eps;
        } else {
            // Два пересечения
            double l = y3p - y3m;
            return (l + (h2 - l) / eps) / h2;
        }
    }

    // Одна точка снаружи, одна внутри: одно пересечение
    double intersectionY = (y3p > min_y && y3p < max_y) ? y3p : y3m;
    double l = (intersectionY == y3p) ? y3p - min_y : max_y - y3m;
    return (l + (h2 - l) / eps) / h2;
}


/**
 * Функция, определяющая коэффициент b[i][j] по положению
 * отрезка относительно кривой(он гарантированно состоит из
 * внутренних узлов сетки).
 * Уравнение эллипса:
 *      x^2 + 4y^2 = 1
 * Внешняя фиктивная область:
 *      (x,y) = [-1.5; 1.5] x [-1; 1]
 * Возвращает:
 *      1 - если отрезок полностью в области D
 *      1 / eps - если отрезок полностью в области D^
 *      l/h2 + (1 - l/h2) / eps - если отрезок пересекает Gamma
 */

double get_bCoeff(const Point& pnt1, const Point& pnt2, const double eps)
{
    bool p1 = inEllips(pnt1.x, pnt1.y), p2 = inEllips(pnt2.x, pnt2.y);

    // Обе точки внутри эллипса
    if (p1 && p2) {
        return 1.;
    } 
    
    // Проверка наличия точек пересечения с эллипсом
    double deg = 1.0 - 4. * sqr(pnt1.y);
    if (deg < 1e-9) { // Обе точки снаружи, без пересечений
        return 1.0 / eps;
    }
    
    double x3p = sqrt(deg), x3m = -sqrt(deg);
    double h1 = std::abs(pnt1.x - pnt2.x);
    double min_x = std::min(pnt1.x, pnt2.x), max_x = std::max(pnt1.x, pnt2.x);
    
    
    // Обе точки снаружи
    if (!p1 && !p2) {
        if (x3p <= min_x || x3m >= max_x) {
            // Без пересечений
            return 1.0 / eps;
        } else {
            // Два пересечения
            double l = x3p - x3m;
            return (l + (h1 - l) / eps) / h1;
        }
    }

    // Одна точка снаружи, одна внутри: одно пересечение
    double intersectionX = (x3p > min_x && x3p < max_x) ? x3p : x3m;
    double l = (intersectionX == x3p) ? x3p - min_x : max_x - x3m;
    return (l + (h1 - l) / eps) / h1;
}





/**
 * Функция, определяющая коэффициент f[i][j] по положению
 * ячейки разбиения П[i][j] относительно кривой(она
 * гарантированно состоит из внутренних узлов сетки).
 * Уравнение эллипса:
 *      x^2 + 4y^2 = 1
 * Внешняя фиктивная область:
 *      (x,y) = [-1.5; 1.5] x [-1; 1]
 * Возвращает:
 *      1 - если ячейка полностью в области D
 *      0 - если ячейка полностью в области D^
 *      S * f(x*, y*) / h1h2 - если ячейка пересекает Gamma
 */
double get_fCoeff(const double x1, const double x2, const double y1, const double y2)
{
    bool p11 = inEllips(x1, y1), p12 = inEllips(x1, y2);
    bool p21 = inEllips(x2, y1), p22 = inEllips(x2, y2);

    // Все точки внутри эллипса
    if (p11 && p12 && p21 && p22) {
        return 1.;
    }

    // Все точки вне эллипса
    if (!p11 && !p12 && !p21 && !p22) {
        return 0.;
    } 

    // Есть пересечение
    double s, h1 = std::abs(x2 - x1), h2 = std::abs(y2 - y1);
    std::vector<Point> vertices;

    if (p11 + p12 + p21 + p22 == 3) // одна точка вне D
    {
        double x4, y4, x5, y5;
        if (!p11 && p12 && p21 && p22) {
            x4 = x1, y4 = sqrt((1. - sqr(x1)) / 4.);
            x5 = -sqrt(1 - 4.* sqr(y1)), y5 = y1;
            vertices.push_back(Point(x2, y1));
            vertices.push_back(Point(x2, y2));
            vertices.push_back(Point(x1, y2));
        } else if (p11 && !p12 && p21 && p22) {
            x4 = x1, y4 = -sqrt((1. - sqr(x1)) / 4.);
            x5 = -sqrt(1 - 4.* sqr(y2)), y5 = y2;
            vertices.push_back(Point(x2, y2));
            vertices.push_back(Point(x2, y1));
            vertices.push_back(Point(x1, y1));
        } else if (p11 && p12 && !p21 && p22) {
            x4 = x2, y4 = sqrt((1. - sqr(x2)) / 4.);
            x5 = sqrt(1 - 4.* sqr(y1)), y5 = y1;
            vertices.push_back(Point(x1, y1));
            vertices.push_back(Point(x1, y2));
            vertices.push_back(Point(x2, y2));
        } else {
            x4 = x2, y4 = -sqrt((1. - sqr(x2)) / 4.);
            x5 = sqrt(1 - 4.* sqr(y2)), y5 = y2;
            vertices.push_back(Point(x1, y2));
            vertices.push_back(Point(x1, y1));
            vertices.push_back(Point(x2, y1));
        }
        vertices.push_back(Point(x4, y4));
        vertices.push_back(Point(x5, y5));    
        return get_Area(vertices) / (h1 * h2);
    }

    if (p11 + p12 + p21 + p22 == 2) // две точки вне D
    {
        double x3, y3, x4, y4;
        if (!p12 && !p22 && p21 && p11) {
            x3 = x2, y3 = -sqrt((1. - sqr(x2)) / 4.);
            x4 = x1, y4 = -sqrt((1. - sqr(x1)) / 4.);
            vertices.push_back(Point(x1, y1));
            vertices.push_back(Point(x2, y1));
        } else if (p12 && !p22 && !p21 && p11) {
            x3 = sqrt(1. - 4. * sqr(y2)), y3 = y2;
            x4 = sqrt(1. - 4. * sqr(y1)), y4 = y1;
            vertices.push_back(Point(x1, y1));
            vertices.push_back(Point(x1, y2));
        } else if (p12 && p22 && !p21 && !p11) {
            x3 = x2, y3 = sqrt((1. - sqr(x2)) / 4.);
            x4 = x1, y4 = sqrt((1. - sqr(x1)) / 4.);
            vertices.push_back(Point(x1, y2));
            vertices.push_back(Point(x2, y2));
        } else {
            x3 = -sqrt(1. - 4. * sqr(y1)), y3 = y1;
            x4 = -sqrt(1. - 4. * sqr(y2)), y4 = y2;
            vertices.push_back(Point(x2, y2));
            vertices.push_back(Point(x2, y1));
        }
        vertices.push_back(Point(x3, y3));   
        vertices.push_back(Point(x4, y4)); 
        return get_Area(vertices) / (h1 * h2);
    }


    // Три точки вне
    double x3, y3, x4, y4;
    if (p11) {
        x3 = sqrt(1. - 4. * sqr(y1)), y3 = y1;
        x4 = x1, y4 = -sqrt((1. - sqr(x1)) / 4.);
        vertices.push_back(Point(x1, y1));
    } else if (p12) {
        x3 = x1, y3 = sqrt((1. - sqr(x1)) / 4.);
        x4 = sqrt(1. - 4. * sqr(y2)), y4 = y2;
        vertices.push_back(Point(x1, y2));
    } else if (p22) {
        x3 = -sqrt(1. - 4. * sqr(y2)), y3 = y2;
        x4 = x2, y4 = sqrt((1. - sqr(x2)) / 4.);
        vertices.push_back(Point(x2, y2));
    } else {
        x3 = x2, y3 = -sqrt((1. - sqr(x2)) / 4.);
        x4 = -sqrt(1. - 4. * sqr(y1)), y4 = y1;
        vertices.push_back(Point(x2, y1));
    }
    vertices.push_back(Point(x3, y3));   
    vertices.push_back(Point(x4, y4)); 
    return get_Area(vertices) / (h1 * h2);
}


