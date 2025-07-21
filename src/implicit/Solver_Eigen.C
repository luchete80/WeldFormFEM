#include "Solver_Eigen.h"
#include "Matrix.h"
#include "Domain_d.h"
#include <iomanip> //setprecision

using namespace std;

namespace MetFEM{
  
void Solver_Eigen::Allocate(){
  m_dof = m_dom->m_node_count * m_dom->m_dim;
  cout << "Allocate for Domain DOFs: "<< m_dof<<endl;
  K.resize(m_dof,m_dof);
  //U.resize(m_dof);
  U = Eigen::VectorXd::Zero(m_dof);
  //R.resize(m_dof);
  R = Eigen::VectorXd::Zero(m_dof);

  // Now K has positive dimensions and you can safely call .sum(), .norm(), etc.
  std::cout << "K.sum() = " << K.sum() << std::endl;
  
}

///// RELIES ON MEMORY SAVED Kelemen.

void Solver_Eigen::assemblyGlobalMatrix() {
    int ndof = m_dom->m_nodxelem * m_dom->m_dim; // DOFs per element
    cout << "ndof "<<ndof <<endl;
    std::vector<T> triplets;
    triplets.reserve(m_dom->m_elem_count * ndof * ndof);

    for (int e = 0; e < m_dom->m_elem_count; ++e) {
        cout << "positioning "<<endl;
        Matrix *Ke = m_dom->m_Kmat[e];
        cout <<"endl"<<endl;
        std::vector<int> global_dofs(ndof);  // <<--- FIXED
        cout <<"set"<<endl;
        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            int node = m_dom->getElemNode(e, a);
            cout << "node "<<node<<endl;
            //int node = m_dom->m_elnod[e*m_dom->m_nodxelem+a];
            for (int i = 0; i < m_dom->m_dim; ++i) {
                global_dofs[a * m_dom->m_dim + i] = node * m_dom->m_dim + i;
            }
        }
        cout << "global dof set"<<endl;
        for (int i = 0; i < ndof; ++i) {  // <<--- FIXED
            int I = global_dofs[i];
            for (int j = 0; j < ndof; ++j) {  // <<--- FIXED
                int J = global_dofs[j];
                cout <<"get val"<<endl;
                double val = Ke->getVal(i, j);
                cout << val<<endl;
                if (val != 0.0) {
                    cout << "I J "<< I <<", "<<J <<", val"<<val<<endl;
                    triplets.emplace_back(I, J, val);
                }
                // Optional: keep or remove this debug print
            }
        }
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
    std::cout << "Matrix mat:\n" << K << std::endl;
    
    //~ cout << "VALUES 9 10 "<<K.coeff(9, 10);
    //~ cout << "VALUES 10 9"<<K.coeff(10,9);
    //~ cout << "VALUES 8 9 "<<K.coeff(8, 9);
    
    
}

//~ void Solver_Eigen::SetBCs(){
  
  //~ for (int dim=0;dim<m_dom->m_dim; dim++){
  //~ par_loop (n,m_dom->bc_count[dim]){
    //~ double val;
    //~ int dof;
    
    //~ if (dim == 0)       {
      //~ dof = m_dom->m_dim*bcx_nod[n]+dim 
      //~ val = m_dom->bcx_val[n];
    //~ } else if (dim == 1)  {m_dom->v[m_dom->m_dim*bcy_nod[n]+dim] = bcy_val[n];
    //~ } else if (dim == 2)  {m_dom->v[m_dom->m_dim*bcz_nod[n]+dim] = bcz_val[n];
    //~ }
  //~ }//node loop
  
//~ }

///// TODO: MODIFY WITH EIGEN MAP FOR PERFORMANCE: 
///// Eigen::Map<Eigen::VectorXi> bc_nodes(m_dom->bcx_nod, m_dom->bc_count[0]);
void Solver_Eigen::applyDirichletBCs() {
    // Phase 1: Identify all BC DOFs (parallel marking)
    std::vector<bool> is_bc_dof(m_dof, false);
    
    #pragma omp parallel for
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        int* nodes = nullptr;
        int count = m_dom->bc_count[dim];
        cout << "DIM "<< dim << " BC Count: "<<count<<endl;
        if (dim == 0) nodes = m_dom->bcx_nod;
        else if (dim == 1) nodes = m_dom->bcy_nod;
        else if (dim == 2) nodes = m_dom->bcz_nod;

        for (int i = 0; i < count; ++i) {
            int dof = nodes[i] * m_dom->m_dim + dim;
            is_bc_dof[dof] = true;  // Thread-safe for distinct DOFs
        }
    }

    // Phase 2: Apply BCs (serial but efficient)
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        int* nodes = nullptr;
        double* values = nullptr;
        int count = m_dom->bc_count[dim];

        if (dim == 0) {
            nodes = m_dom->bcx_nod;
            values = m_dom->bcx_val;
        } else if (dim == 1) {
            nodes = m_dom->bcy_nod;
            values = m_dom->bcy_val;
        } else if (dim == 2) {
            nodes = m_dom->bcz_nod;
            values = m_dom->bcz_val;
        }

        for (int i = 0; i < count; ++i) {
            int dof = nodes[i] * m_dom->m_dim + dim;
            double value = values[i];
            cout << "dof: "<<dof<<", "<<value<<endl;
            
       // Clear row
       //// THE SAME AS FOR COLUMN IMPLEMENTED WITH ROW DOES NOT WORK.
        for (int k = K.outerIndexPtr()[dof]; k < K.outerIndexPtr()[dof + 1]; ++k) {
            if (K.innerIndexPtr()[k] != dof)
                K.valuePtr()[k] = 0.0;
        }

        // Clear column
        for (int col = 0; col < K.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it) {
                if (it.row() == dof && it.col() != dof) {
                    cout << "clear row "<<it.row() <<", "<<it.col()<<endl; 
                    it.valueRef() = 0.0;
                }
            }
        }

            // Set diagonal and RHS
            K.coeffRef(dof, dof) = 1.0;
            R[dof] = value;
        }
    }
    
    
    //~ R[dof] = value;
    //~ R[dof] = value;
    //~ R[dof] = value;
    
    // Optional but recommended for performance
    K.makeCompressed();
    std::cout << "Matrix mat:\n" << K << std::endl;
    
}

void Solver_Eigen::SetRDOF(const int &dof, const double &val){
  R[dof] = val;
  
}


int Solver_Eigen::Solve(){


    std::cout << "Matrix mat:\n" << K << std::endl;
    std::cout << "Vector mat:\n" << R << std::endl;
    
    //R << 1, -2, 0;
    std::cout << "K norm: " << K.norm() << std::endl;
    cout << "Analyzing pattern"<<endl;
    solver.analyzePattern(K);
    solver.setPivotThreshold(1e-6);  // Allow smaller pivots

    cout << "Done. "<<endl;
    cout << "Factorizing"<<endl;
    solver.factorize(K);
    

    if(solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed\n";
        return -1;
    }

    U = solver.solve(R);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Solving failed\n";
        return -1;
    }

    if (solver.info() != Eigen::Success) {
      std::cerr << "Factorization failed, info code: " << solver.info() << std::endl;
      return -1;
    }

    std::cout << "Solution:\n" << std::setprecision(8) <<U << std::endl; // Should print [1, -1, -2]

    return 0;
  
  
}

    void Solver_Eigen::assembleElement(int e, /*const*/ Matrix& Ke) {
        std::vector<int> global_dofs(m_dom->m_nodxelem * m_dom->m_dim);
        
        // Compute global DOFs for this element
        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            const int node = m_dom->getElemNode(e, a);
            for (int i = 0; i < m_dom->m_dim; ++i) {
                global_dofs[a * m_dom->m_dim + i] = node * m_dom->m_dim + i;
            }
        }
        
        addElementToTriplets(Ke, global_dofs, m_triplets);
    }

    // Finalize the assembly
    void Solver_Eigen::finalizeAssembly() {
        K.setFromTriplets(m_triplets.begin(), m_triplets.end());
        K.makeCompressed();
        m_triplets.clear(); // Free memory
    }

};
