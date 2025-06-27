#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>

int main() {
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::Triplet<double> T;

    int n = 3;
    std::vector<T> triplets;
    triplets.reserve(7);

    // Build the matrix A (3x3)
    triplets.push_back(T(0,0,3)); triplets.push_back(T(0,1,2)); triplets.push_back(T(0,2,-1));
    triplets.push_back(T(1,0,2)); triplets.push_back(T(1,1,-2)); triplets.push_back(T(1,2,4));
    triplets.push_back(T(2,0,-1)); triplets.push_back(T(2,1,0.5)); triplets.push_back(T(2,2,-1));

    SpMat A(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd b(n);
    b << 1, -2, 0;

    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    if(solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed\n";
        return -1;
    }

    Eigen::VectorXd x = solver.solve(b);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Solving failed\n";
        return -1;
    }

    std::cout << "Solution:\n" << x << std::endl; // Should print [1, -1, -2]

    return 0;
}
