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
  
  // THIS REQUIRES ALL ELEMENT MATRICES TO BE STORED!!!
  virtual void assemblyGlobalMatrix();
  ////////////////////////////////////
  virtual void Allocate();
  virtual void applyDirichletBCs();
  virtual void SetRDOF(const int &, const double &val);

    // Main assembly function that can be called element-by-element
    void beginAssembly() {
        const int ndof = m_dom->m_nodxelem * m_dom->m_dim;
        const size_t estimated_nnz = m_dom->m_elem_count * ndof * ndof;
        m_triplets.clear();
        m_triplets.reserve(estimated_nnz);
    }

    // Add a single element to the assembly (NOT NEED TO BE STORED!)
    void assembleElement(int e, const Matrix& Ke);

    // Finalize the assembly
    void finalizeAssembly();
    /////// USAGE:
    // In your time step loop:
    // // solver.beginAssembly();

    // // for (int e = 0; e < m_dom->m_elem_count; ++e) {
        // // // 1. Compute your element stiffness matrix (Kmat + Kgeo)
        // // Matrix Ke = computeElementStiffness(e); 
        
        // // // 2. Apply BCs to this element matrix if needed
        // // applyBCsToElement(e, Ke);
        
        // // // 3. Add mass scaling if needed
        // // addMassScaling(e, Ke);
        
        // // // 4. Add to global assembly
        // // solver.assembleElement(e, Ke);
    // // }

    // // solver.finalizeAssembly();
    
    
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
  std::vector<T> m_triplets;  // Almacenamiento temporal para los triplets durante el ensamblaje



private:
    // Helper function to add element matrix contributions to triplets
    inline void addElementToTriplets(/*const */Matrix& Ke, 
                                   const std::vector<int>& global_dofs,
                                   std::vector<T>& triplets) {
        const int ndof = global_dofs.size();
        
        for (int i = 0; i < ndof; ++i) {
            const int I = global_dofs[i];
            for (int j = 0; j < ndof; ++j) {
                const double val = Ke.getVal(i, j);
                if (val != 0.0) {
                    triplets.emplace_back(I, global_dofs[j], val);
                }
            }
        }
    }
  
};



};
