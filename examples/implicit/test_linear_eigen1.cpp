#include <iostream>
#include <Eigen/Dense>

int main() {
    Eigen::Matrix3d A;
    A << 3, 2, -1,
         2, -2, 4,
        -1, 0.5, -1;

    Eigen::Vector3d b(1, -2, 0);

    // Solve Ax = b
    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);

    std::cout << "Solution x:\n" << x << std::endl;

    // For the second system:
    Eigen::Matrix2d A2;
    A2 << 2, 1,
          1, 3;

    Eigen::Vector2d b2(1, 2);

    Eigen::Vector2d x2 = A2.colPivHouseholderQr().solve(b2);

    std::cout << "Solution x2:\n" << x2 << std::endl;

    return 0;
}
