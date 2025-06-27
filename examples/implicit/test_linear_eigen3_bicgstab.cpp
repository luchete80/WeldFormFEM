#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

int main() {
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::Triplet<double> T;

    int n = 3;
    std::vector<T> triplets;

    // Fill matrix A (example matrix, can be non-symmetric)
    triplets.push_back(T(0, 0, 3));
    triplets.push_back(T(0, 1, 2));
    triplets.push_back(T(0, 2, -1));
    triplets.push_back(T(1, 0, 2));
    triplets.push_back(T(1, 1, -2));
    triplets.push_back(T(1, 2, 4));
    triplets.push_back(T(2, 0, -1));
    triplets.push_back(T(2, 1, 0.5));
    triplets.push_back(T(2, 2, -1));

    SpMat A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd b(n);
    b << 1, -2, 0;

    // Create BiCGSTAB solver and compute the factorization
    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(A);

    if(solver.info() != Eigen::Success) {
        std::cout << "Decomposition failed" << std::endl;
        return -1;
    }

    // Solve Ax = b
    Eigen::VectorXd x = solver.solve(b);

    if(solver.info() != Eigen::Success) {
        std::cout << "Solver failed" << std::endl;
        return -1;
    }

    std::cout << "Solution x:\n" << x << std::endl;

    return 0;
}
