#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "util.h"

namespace backend {

struct Face {
    int v0;
    int v1;	
    Face() {}
    Face(int v0, int v1) :v0(v0), v1(v1) {}
};

// n1, n2: cell numbers
MatrixXd computeSDF_cuboid(const int& n1, const int& n2, const Vector2i& left_corner, const Vector2i& cuboid_size);
double computeSDF_cuboid(const Vector2d& pos, const Vector2d& cuboid_size);    
MatrixXd computeSDF_circle(const int& n1, const int& n2, const Vector2d& center, const double& r);

}

#endif