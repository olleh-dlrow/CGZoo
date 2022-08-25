/*
c++ eigen_test.cpp -std=c++11 -o eigen_test && ./eigen_test
*/
#include<iostream>
#include"Eigen/Core"
#include"Eigen/Dense"
using namespace Eigen;

int main() {
    Vector3f v3(4, 5, 6);
    std::cout << v3 * Vector3f(1, 2, 3);
    return 0;
}