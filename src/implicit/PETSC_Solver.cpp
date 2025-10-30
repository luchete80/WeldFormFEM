

/************************************************************************

Copyright 2012-2013 Luciano Buglioni

Contact:
    luciano.buglioni@gmail.com

This file is a part of FluxSol

For a copy of the GNU General Public License,
see <http://www.gnu.org/licenses/>.

*************************************************************************/

#include "PETSC_Solver.h"

namespace MetFEM
{
  
  ////////////////////////// ALLOCATE
  
  std::vector<int> nnz_per_row(num_dofs, 0);

// Loop over elements
for (int e = 0; e < m_elem_count; ++e) {
    std::vector<int> global_dofs(m_nodxelem * m_dim);
    for (int a = 0; a < m_nodxelem; ++a) {
        int node = getElemNode(e, a);
        for (int i = 0; i < m_dim; ++i) {
            global_dofs[a * m_dim + i] = node * m_dim + i;
        }
    }

    for (int i = 0; i < global_dofs.size(); ++i) {
        int row = global_dofs[i];
        for (int j = 0; j < global_dofs.size(); ++j) {
            int col = global_dofs[j];
            if (col != row) nnz_per_row[row]++;
        }
        nnz_per_row[row]++; // Diagonal
    }
}
  
  
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, num_dofs, num_dofs);
  MatSetFromOptions(A);
  MatSeqAIJSetPreallocation(A, 0, nnz_per_row.data());
  MatMPIAIJSetPreallocation(A, 0, nnz_per_row.data(), 0, NULL); // optional if MP



    ///////////////////////////////////////////////
  /////////////////////////////// SET VALUES
  
  
    for (int e = 0; e < m_elem_count; ++e) {
      const Matrix& Ke = *(m_Kmat[e]);

      std::vector<int> global_dofs(m_nodxelem * m_dim);
      for (int a = 0; a < m_nodxelem; ++a) {
          int node = getElemNode(e, a);
          for (int i = 0; i < m_dim; ++i) {
              global_dofs[a * m_dim + i] = node * m_dim + i;
          }
      }

      std::vector<PetscScalar> values(global_dofs.size() * global_dofs.size());

      for (int i = 0; i < global_dofs.size(); ++i) {
          for (int j = 0; j < global_dofs.size(); ++j) {
              values[i * global_dofs.size() + j] = Ke(i, j);
          }
      }

      MatSetValues(A,
          global_dofs.size(), global_dofs.data(),  // Rows
          global_dofs.size(), global_dofs.data(),  // Cols
          values.data(), ADD_VALUES);
  }

  
  
  
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  
  
  ////ASSEMBLY
  
void Solver_Eigen::assemblyGlobalMatrix() {
  
  for (int e = 0; e < m_elem_count; ++e) {
    const Matrix& Ke = *(m_Kmat[e]);

    std::vector<int> global_dofs(m_nodxelem * m_dim);
    for (int a = 0; a < m_nodxelem; ++a) {
        int node = getElemNode(e, a);
        for (int i = 0; i < m_dim; ++i) {
            global_dofs[a * m_dim + i] = node * m_dim + i;
        }
    }

    std::vector<PetscScalar> values(global_dofs.size() * global_dofs.size());

    for (int i = 0; i < global_dofs.size(); ++i) {
        for (int j = 0; j < global_dofs.size(); ++j) {
            values[i * global_dofs.size() + j] = Ke(i, j);
        }
    }

    MatSetValues(A,
        global_dofs.size(), global_dofs.data(),  // Rows
        global_dofs.size(), global_dofs.data(),  // Cols
        values.data(), ADD_VALUES);
  }

}


} //METFEM

