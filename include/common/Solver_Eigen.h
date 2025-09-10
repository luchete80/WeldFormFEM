#include "Solver.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>

#include <unordered_map>  // For triplet storage

#include "defs.h"
#include "Domain_d.h"

class Matrix;

using namespace std;

namespace MetFEM{
  
class Solver_Eigen:
public Solver{

public: 
  virtual int Solve();
  Solver_Eigen(){
    incremental_bcs = std::vector<std::pair<int, double>>(); // Construcción explícita
    incremental_bcs.reserve(100); // Reservar memoria inicial
    }
  virtual ~Solver_Eigen(){

    }

  inline const double & getU(int node, int dim) const {
      return U[node * m_dom->m_dim + dim]; 
  }

  inline const double & getR(int node, int dim) const {
      return R[node * m_dom->m_dim + dim]; 
  }
  
  inline void addToU(int node, int dim, double delta) {
      m_dom->u[node * m_dom->m_dim + dim] += delta;
      m_dom->x[node * m_dom->m_dim + dim] += delta; // Update position as well
  }
  
  void printK(){
    std::cout << "Vector mat:\n" << K << std::endl;
  }

  void printR(){
    std::cout << "Vector mat:\n" << R << std::endl;
  }
  
  
  void setZero(){
    R.setZero();
    K.setZero(); //May be not necesary since is made from triplets
    }
  // THIS REQUIRES ALL ELEMENT MATRICES TO BE STORED!!!
  virtual void assemblyGlobalMatrix();
  ////////////////////////////////////
  virtual void Allocate();
  virtual void applyDirichletBCs();
  virtual void SetRDOF(const int &, const double &val);

  // void getElementDisplacements(int e, Eigen::VectorXd& u_e) {
      // u_e.resize(m_dom->m_nodxelem * m_dom->m_dim);
      // for (int a = 0; a < m_dom->m_nodxelem; ++a) {
          // int node = m_dom->getElemNode(e, a);
          // for (int i = 0; i < m_dom->m_dim; ++i) {
              // u_e(a * m_dom->m_dim + i) = U(node * m_dom->m_dim + i);
          // }
      // }
  // }
  //inline double getUNode(const int&n,){return U[m_dim*]};

    // Main assembly function that can be called element-by-element
    void beginAssembly() {
        const int ndof = m_dom->m_nodxelem * m_dom->m_dim;
        const size_t estimated_nnz = m_dom->m_elem_count * ndof * ndof;
        m_triplets.clear();
        m_triplets.reserve(estimated_nnz);
    }

    // Add a single element to the assembly (NOT NEED TO BE STORED!)
    void assembleElement(int e, /*const*/ Matrix& Ke);

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
    virtual double getRNorm(){return R.norm();}
    virtual double getUNorm(){return U.norm();}
    ///virtual double getCondNumber(){return K.norm() * K.inverse().norm();} //PROHIBITIVE
    void assembleResidual(int e, /*const */ Matrix& Re) {
        const int ndof = m_dom->m_nodxelem * m_dom->m_dim;
        std::vector<int> global_dofs(ndof);
        
        // Precompute DOF mapping - helps compiler optimize
        //#pragma unroll(4) // Encourage loop unrolling for small elements
        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            const int node = m_dom->getElemNode(e, a);
            const int base = a * m_dom->m_dim;
            for (int i = 0; i < m_dom->m_dim; ++i) {
                global_dofs[base + i] = node * m_dom->m_dim + i;
            }
        }

        // Atomic adds for thread safety in parallel loops
        //#pragma omp parallel for if(ndof > 128) // Threshold for parallelization
        for (int i = 0; i < ndof; ++i) {
            #pragma omp atomic update
            R(global_dofs[i]) += Re.getVal(i, 0);
        }
  }

  inline void addToR(int global_dof, double value) {
      #pragma omp atomic update  // For thread safety if using OpenMP
      R(global_dof) += value;
  }
  virtual void assembleContactStiffness(double kn, double dt);
  
  virtual void setDirichletBC(int dof, double delta_value);
  void applyIncrementalBCs();

protected:
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Triplet<double> T;

  
  SpMat K;
  Eigen::VectorXd R;
  
  Eigen::VectorXd U;
  
  Eigen::SparseLU<SpMat> solver;
  std::vector<T> m_triplets;  // Almacenamiento temporal para los triplets durante el ensamblaje
  
  std::vector<std::pair<int, double>> incremental_bcs;


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
