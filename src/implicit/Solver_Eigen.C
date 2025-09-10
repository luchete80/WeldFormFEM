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
    //cout << "ndof "<<ndof <<endl;
    std::vector<T> triplets;
    triplets.reserve(m_dom->m_elem_count * ndof * ndof);

    for (int e = 0; e < m_dom->m_elem_count; ++e) {
        //cout << "positioning "<<endl;
        Matrix *Ke = m_dom->m_Kmat[e];
        //cout <<"endl"<<endl;
        std::vector<int> global_dofs(ndof);  // <<--- FIXED
        //cout <<"set"<<endl;
        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            int node = m_dom->getElemNode(e, a);
            //cout << "node "<<node<<endl;
            //int node = m_dom->m_elnod[e*m_dom->m_nodxelem+a];
            for (int i = 0; i < m_dom->m_dim; ++i) {
                global_dofs[a * m_dom->m_dim + i] = node * m_dom->m_dim + i;
            }
        }
        //cout << "global dof set"<<endl;
        for (int i = 0; i < ndof; ++i) {  // <<--- FIXED
            int I = global_dofs[i];
            for (int j = 0; j < ndof; ++j) {  // <<--- FIXED
                int J = global_dofs[j];
                //cout <<"get val"<<endl;
                double val = Ke->getVal(i, j);
                cout << val<<endl;
                if (val != 0.0) {
                    //cout << "I J "<< I <<", "<<J <<", val"<<val<<endl;
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
//// ERROR: THIS DOES NOT DEAL WITH NONZERO VALS
///// Eigen::Map<Eigen::VectorXi> bc_nodes(m_dom->bcx_nod, m_dom->bc_count[0]);
// void Solver_Eigen::applyDirichletBCs() {
    // // Phase 1: Identify all BC DOFs (parallel marking)
    // std::vector<bool> is_bc_dof(m_dof, false);
    
    // #pragma omp parallel for
    // for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        // int* nodes = nullptr;
        // int count = m_dom->bc_count[dim];
        // cout << "DIM "<< dim << " BC Count: "<<count<<endl;
        // if (dim == 0) nodes = m_dom->bcx_nod;
        // else if (dim == 1) nodes = m_dom->bcy_nod;
        // else if (dim == 2) nodes = m_dom->bcz_nod;

        // for (int i = 0; i < count; ++i) {
            // int dof = nodes[i] * m_dom->m_dim + dim;
            // is_bc_dof[dof] = true;  // Thread-safe for distinct DOFs
        // }
    // }

    // // Phase 2: Apply BCs (serial but efficient)
    // for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        // int* nodes = nullptr;
        // double* values = nullptr;
        // int count = m_dom->bc_count[dim];

        // if (dim == 0) {
            // nodes = m_dom->bcx_nod;
            // values = m_dom->bcx_val;
        // } else if (dim == 1) {
            // nodes = m_dom->bcy_nod;
            // values = m_dom->bcy_val;
        // } else if (dim == 2) {
            // nodes = m_dom->bcz_nod;
            // values = m_dom->bcz_val;
        // }

        // for (int i = 0; i < count; ++i) {
            // int dof = nodes[i] * m_dom->m_dim + dim;
            // double value = values[i];
            // cout << "dof: "<<dof<<", "<<value<<endl;
            
       // // Clear row
       // //// THE SAME AS FOR COLUMN IMPLEMENTED WITH ROW DOES NOT WORK.
        // for (int k = K.outerIndexPtr()[dof]; k < K.outerIndexPtr()[dof + 1]; ++k) {
            // if (K.innerIndexPtr()[k] != dof)
                // K.valuePtr()[k] = 0.0;
        // }

        // // Clear column
        // for (int col = 0; col < K.outerSize(); ++col) {
            // for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it) {
                // if (it.row() == dof && it.col() != dof) {
                    // //cout << "clear row "<<it.row() <<", "<<it.col()<<endl; 
                    // it.valueRef() = 0.0;
                // }
            // }
        // }

            // // Set diagonal and RHS
            // K.coeffRef(dof, dof) = 1.0;
            // R[dof] = value;
        // }
    // }
    
    
    // //~ R[dof] = value;
    // //~ R[dof] = value;
    // //~ R[dof] = value;
    
    // // Optional but recommended for performance
    // K.makeCompressed();
    // std::cout << "Matrix mat:\n" << K << std::endl;
    
// }

void Solver_Eigen::applyDirichletBCs() {
    // Phase 1: Mark all BC DOFs
    std::vector<std::pair<int, double>> bc_dofs;
    
    // Pre-count total BCs
    int total_bcs = 0;
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        total_bcs += m_dom->bc_count[dim];
    }
    bc_dofs.reserve(total_bcs);

    // Collect all BC DOFs with their values
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        const int count = m_dom->bc_count[dim];
        const int* nodes = nullptr;
        const double* values = nullptr;

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
            const int dof = nodes[i] * m_dom->m_dim + dim;
            bc_dofs.emplace_back(dof, values[i]);
            cout <<"BC value "<<i << ": "<<values[i]<<endl;
        }
    }

    // Create a mapping for quick lookup
    std::unordered_map<int, double> bc_map(bc_dofs.begin(), bc_dofs.end());

    // Phase 2: Adjust RHS first (before modifying K)
    // This is more efficient and avoids issues with modified matrix
    for (const auto& [dof, value] : bc_dofs) {
        // For each BC DOF, subtract K_{i,dof} * value from RHS[i] for all i != dof
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, dof); it; ++it) {
            if (it.row() != dof) {
                R[it.row()] -= it.value() * value;
            }
        }
    }

    // Phase 3: Modify matrix K
    K.makeCompressed();
    
    for (const auto& [dof, value] : bc_dofs) {
        // Clear row (set off-diagonals to zero)
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, dof); it; ++it) {
            if (it.col() != dof) {
                it.valueRef() = 0.0;
            } else {
                it.valueRef() = 1.0; // Set diagonal to 1
            }
        }

        // Clear column (set off-diagonals to zero)
        // This is tricky with Eigen's column-major storage
        // Better approach: use a different strategy
        for (int col = 0; col < K.outerSize(); ++col) {
            if (col == dof) continue; // Skip the diagonal column
            
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it) {
                if (it.row() == dof) {
                    it.valueRef() = 0.0;
                }
            }
        }

        // Set RHS for this DOF
        R[dof] = value;
    }
    
    std::cout << "R after  BCs:\n" << R << std::endl;

    K.makeCompressed(); // Recompress after modifications
}

void Solver_Eigen::SetRDOF(const int &dof, const double &val){
  R[dof] = val;
  
}


int Solver_Eigen::Solve(){


    // std::cout << "Matrix mat:\n" << K << std::endl;
    // std::cout << "Vector mat:\n" << R << std::endl;
    
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

    //std::cout << "Solution:\n" << std::setprecision(8) <<U << std::endl; // Should print [1, -1, -2]

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


// SIMPLEST CASE, NODE TO RIGID
// 

void Solver_Eigen::assembleContactStiffness(double kn, double dt) {
    for (int i = 0; i < m_dom->m_node_count; ++i) {
        if (m_dom->m_mesh_in_contact[i]>-1) continue;  // marcás qué nodos están en contacto

        // normal en el nodo (ya la calculaste en contforce)
        //Eigen::Vector3d n = m_dom->contactNormal[i]; 
        Eigen::Vector3d n (m_dom->contforce[3*i], m_dom->contforce[3*i+1],m_dom->contforce[3*i+2]);
        n /= n.norm();
        Eigen::Matrix3d nnT = n * n.transpose();   // (n ⊗ n)

        // bloque de rigidez de contacto
        Eigen::Matrix3d Kc = (kn * dt) * nnT;

        // ensamblar en la matriz global
        int base = i * m_dom->m_dim;  // desplazamiento DOF del nodo i
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                m_triplets.emplace_back(base + r, base + c, Kc(r,c));
            }
        }
    }
}


    // Finalize the assembly
    void Solver_Eigen::finalizeAssembly() {
        K.setFromTriplets(m_triplets.begin(), m_triplets.end());
        K.makeCompressed();
        m_triplets.clear(); // Free memory
    }
    
    void Solver_Eigen::setDirichletBC(int dof, double delta_value) {
    incremental_bcs.emplace_back(dof, delta_value);
    cout << "adding vector "<<endl;
    //incremental_bcs.push_back(std::make_pair(dof, delta_value));
}


void Solver_Eigen::applyIncrementalBCs() {
    // Aplicar solo BCs incrementales (sin las fijas tradicionales)
    if (incremental_bcs.empty()) return;
    
    // 1. Ajustar RHS
    for (const auto& [dof, delta_value] : incremental_bcs) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, dof); it; ++it) {
            if (it.row() != dof) {
                R[it.row()] -= it.value() * delta_value;
            }
        }
    }

    // 2. Modificar matriz K
    K.makeCompressed();
    
    for (const auto& [dof, delta_value] : incremental_bcs) {
        // Limpiar fila y columna
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, dof); it; ++it) {
            if (it.col() != dof) {
                it.valueRef() = 0.0;
            } else {
                it.valueRef() = 1.0;
            }
        }

        for (int col = 0; col < K.outerSize(); ++col) {
            if (col == dof) continue;
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it) {
                if (it.row() == dof) {
                    it.valueRef() = 0.0;
                }
            }
        }

        R[dof] = delta_value;
    }
    
    incremental_bcs.clear();
    K.makeCompressed();
}

};
