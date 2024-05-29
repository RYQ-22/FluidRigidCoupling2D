#ifndef CONFIG_H
#define CONFIG_H

// Commonly used std headers.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>

// Eigen headers.
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "unsupported/Eigen/Polynomials"
#include "unsupported/Eigen/MatrixFunctions"

// OpenCV
#include <opencv2/opencv.hpp>

namespace backend {

// Define constant
const double PI = 3.1415926535897932384626433832795;
const double G = 9.8;

// Alias
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using MatrixXi = Eigen::MatrixXi;
using Vector3d = Eigen::Vector3d;
using Vector2d = Eigen::Vector2d;
using Matrix2d = Eigen::Matrix2d;
using Matrix3d = Eigen::Matrix3d;
using Vector3i = Eigen::Vector3i;
using Vector2i = Eigen::Vector2i;

using FuncPtr = double (*)(const Vector2d&);

// Path
const std::string projectPath = std::string("/home/ricky/Documents/MyProject/FluidRigidCoupling2D");
const std::string pythonPath = std::string("/home/ricky/anaconda3/envs/myenv/bin/python");

}

#endif