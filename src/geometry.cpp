#include "geometry.h"

namespace backend {

MatrixXd computeSDF_cuboid(const int& n1, const int& n2, const Vector2i& left_corner, const Vector2i& cuboid_size) {
    /*    
        r4 --- r3
        |       |
        |       |
        r1 --- r2
    */
    Vector2i r1, r2, r3, r4;
    r1 = left_corner;
    r2 = left_corner + Vector2i(1, 0) * cuboid_size(0);
    r3 = left_corner + cuboid_size;
    r4 = left_corner + Vector2i(0, 1) * cuboid_size(1);    

    MatrixXd phi(n1+1, n2+1);
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) {
        if (i <= r1(0) && j <= r1(1)) {            
            phi(i, j) = ((r1 - Vector2i(i, j)).cast<double>()).norm();            
        }
        else if (i >= r2(0) && j <= r2(1)) {
            phi(i, j) = ((r2 - Vector2i(i, j)).cast<double>()).norm();            
        }
        else if (i >= r3(0) && j >= r3(1)) {
            phi(i, j) = ((r3 - Vector2i(i, j)).cast<double>()).norm();            
        }
        else if (i <= r4(0) && j >= r4(1)) {
            phi(i, j) = ((r4 - Vector2i(i, j)).cast<double>()).norm();            
        }        
        else if (i < r1(0) && j > r1(1) && j < r4(1)) {
            phi(i, j) = static_cast<double>(r1(0) - i);
        }
        else if (i > r2(0) && j > r2(1) && j < r3(1)) {
            phi(i, j) = static_cast<double>(i - r2(0));
        }
        else if (j < r1(1) && i > r1(0) && i < r2(0)) {
            phi(i, j) = static_cast<double>(r1(1) - j);
        }
        else if (j > r4(1) && i > r4(0) && i < r3(0)) {
            phi(i, j) = static_cast<double>(j - r4(1));
        }
        else {
            phi(i, j) = -static_cast<double>(min(i-r1(0), r2(0)-i, j-r1(1), r4(1)-j));
        }
    }
    return phi;
}

double computeSDF_cuboid(const Vector2d& pos, const Vector2d& cuboid_size) {
    /*    
        r4 --- r3
        |       |
        |       |
        r1 --- r2
    */
    Vector2d r1, r2, r3, r4;
    r1 = Vector2d(-cuboid_size(0) / 2, -cuboid_size(1) / 2);
    r2 = Vector2d(cuboid_size(0) / 2, -cuboid_size(1) / 2);
    r3 = Vector2d(cuboid_size(0) / 2, cuboid_size(1) / 2);
    r4 = Vector2d(-cuboid_size(0) / 2, cuboid_size(1) / 2);    
    double i = pos(0), j = pos(1);
    double phi;
    if (i <= r1(0) && j <= r1(1)) {            
        phi = (r1 - pos).norm();            
    }
    else if (i >= r2(0) && j <= r2(1)) {
        phi = (r2 - pos).norm();            
    }
    else if (i >= r3(0) && j >= r3(1)) {
        phi = (r3 - pos).norm();            
    }
    else if (i <= r4(0) && j >= r4(1)) {
        phi = (r4 - pos).norm();            
    }        
    else if (i < r1(0) && j > r1(1) && j < r4(1)) {
        phi = r1(0) - i;
    }
    else if (i > r2(0) && j > r2(1) && j < r3(1)) {
        phi = i - r2(0);
    }
    else if (j < r1(1) && i > r1(0) && i < r2(0)) {
        phi = r1(1) - j;
    }
    else if (j > r4(1) && i > r4(0) && i < r3(0)) {
        phi = j - r4(1);
    }
    else {
        phi = -min(i-r1(0), r2(0)-i, j-r1(1), r4(1)-j);
    }
    return phi;
}

MatrixXd computeSDF_circle(const int& n1, const int& n2, const Vector2d& center, const double& r) {
    MatrixXd phi(n1+1, n2+1);
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) {
        Vector2d pos(i, j);
        phi(i, j) = (pos - center).norm() - r;
    }
    return phi;
}

}