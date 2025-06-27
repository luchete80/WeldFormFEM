#include "Solver.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>

#include <unordered_map>  // For triplet storage

#include "defs.h"

using namespace std;

namespace MetFEM{
  
class Solver_Eigen:
public Solver{

public: 
  virtual int Solve();
  virtual ~Solver_Eigen(){}
  virtual void assemblyGlobalMatrix();
  virtual void Allocate();
  virtual void applyDirichletBCs();
  virtual void SetRDOF(const int &, const double &val);

    //~ void setSolverDirect() { 
        //~ solver.reset(new Eigen::SparseLU<SpMat>());
    //~ }
    
    //~ void setSolverIterative(double tol = 1e-8) {
        //~ auto iterative = new Eigen::ConjugateGradient<SpMat>();
        //~ iterative->setTolerance(tol);
        //~ solver = std::move(iterative);
    //~ }
    
protected:
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Triplet<double> T;

  
  SpMat K;
  Eigen::VectorXd R;
  
  Eigen::VectorXd U;
  
  Eigen::SparseLU<SpMat> solver;
  
};



};
