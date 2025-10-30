/*************************************************************************/
/*  Solver.C                                                     */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/


#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"

#include "Matrix.h"

#include "Solver_Eigen.h"


//~ for each element e:
  //~ // 1. Calcute deformation gradient F
  //~ F = Σ [x_a ⊗ ∇N_a]  // Current config
  //~ F_old = Σ [x_old_a ⊗ ∇N_a]  // Prefvious config
  
  //~ // 2. Calcular tensor de velocidad de deformación ε
  //~ L = (1/Δt) * (F·F_old⁻¹ - I)  // Gradiente de velocidad
  //~ ε = 0.5*(L + Lᵀ)  // Parte simétrica (tasa de deformación)
  
  //~ // 3. Calcular tensiones
  //~ σ = D:ε  // Ley elástica
  
  //~ // 4. Ensamblar fuerzas internas
  //~ f_int = V_e * Bᵀ · σ  // B: matriz deformación-desplazamiento
  
  //~ // 5. Ensamblar matriz de rigidez
  //~ K_mat = V_e * Bᵀ · D · B  // Parte material
  //~ K_geo = Parte geométrica (σ-dependiente)
  
  //~ // 6. Aplicar condiciones de contorno
  //~ // 7. Resolver sistema local: K · Δu = f_int
  //~ // 8. Actualizar desplazamientos
  
  
  
  
// #ifdef BUILD_REMESH
// #include "ReMesher.h"
// #endif

// #include <Eigen/Dense>
// Matrix Matrix::SolveEigen(const Matrix& b) const {
    // Eigen::MatrixXd A_eigen(rows, cols);
    // Eigen::VectorXd b_eigen(rows);
    
    // // Copiar datos a Eigen
    // for (int i = 0; i < rows; ++i) {
        // b_eigen(i) = b.getVal(i, 0);
        // for (int j = 0; j < cols; ++j) {
            // A_eigen(i, j) = getVal(i, j);
        // }
    // }
    
    // // Resolver con Cholesky (LLT)
    // Eigen::VectorXd x_eigen = A_eigen.llt().solve(b_eigen);
    
    // // Devolver resultado
    // Matrix x(rows, 1);
    // for (int i = 0; i < rows; ++i) {
        // x.Set(i, 0, x_eigen(i));
    // }
    // return x;
// }

/// FOR PRODUCTION (SLOW COMPILATION)
// // Ejemplo con Eigen (más rápido y robusto)
// #include <Eigen/Dense>

// //////// PUT CONST
// Matrix SolveWithEigen( Matrix& A,  Matrix& b) {
    // Eigen::MatrixXd Aeig(A.m_row, A.m_col);
    // Eigen::VectorXd beig(b.m_row);
    
    // // Copiar datos
    // for (int i = 0; i < A.m_row; ++i) {
        // beig(i) = b.getVal(i, 0);
        // for (int j = 0; j < A.m_col; ++j) {
            // Aeig(i, j) = A(i, j);
        // }
    // }
    
    // // Resolver
    // Eigen::VectorXd xeig = Aeig.partialPivLu().solve(beig);
    
    // Matrix x(A.m_row, 1);
    // for (int i = 0; i < A.m_row; ++i) {
        // x(i, 0) = xeig(i);
    // }
    // return x;
// }

using namespace std;

// void SolveCholesky(Matrix& K, double* rhs, double* result) {
    // int n = K.m_row;
    // Matrix L(n, n);

    // // Cholesky decomposition: K = L * Lᵗ
    // for (int i = 0; i < n; ++i) {
        // for (int j = 0; j <= i; ++j) {
            // double sum = K.getVal(i, j);
            // for (int k = 0; k < j; ++k)
                // sum -= L.getVal(i, k) * L.getVal(j, k);

            // if (i == j)
                // L.Set(i, j, sqrt(sum));
            // else
                // L.Set(i, j, sum / L.getVal(j, j));
        // }
    // }

    // // Forward substitution: L * y = rhs
    // std::vector<double> y(n);
    // for (int i = 0; i < n; ++i) {
        // double sum = rhs[i];
        // for (int j = 0; j < i; ++j)
            // sum -= L.getVal(i, j) * y[j];
        // y[i] = sum / L.getVal(i, i);
    // }

    // // Backward substitution: Lᵗ * x = y
    // for (int i = n - 1; i >= 0; --i) {
        // double sum = y[i];
        // for (int j = i + 1; j < n; ++j)
            // sum -= L.getVal(j, i) * result[j];
        // result[i] = sum / L.getVal(i, i);
    // }
// }

namespace MetFEM{
  

////////////////////////////////////////
/////// ONLY FOR TESTING  ELASTIC
////////////////////////////////////////
void host_ Domain_d::ElasticIncSolve(){
  WallTimer timer;

  AssignMatAddress();

  cout << "done"<<endl;


  cout << "Imposing BCS"<<endl;
  
  
  // InitValues();
  
  // for (int d=0;d<m_dim;d++){
    
      // for (int n=0;n<m_node_count*m_dim;n++){
        // v[n]=a[n]=u[n]=0.0;
      // }
       // ImposeBCV(d);
   
  // }

  calcElemJAndDerivatives();
  
  //FOR PRESSURE ANP
  //CalcElemInitialVol(); //ALSO CALC VOL

  CalcElemVol();

  //IMPLICIT 
  CalcMaterialStiffElementMatrix();


  //printf("calc dens\n");
  //calcElemDensity();

  //CalcNodalVol(); //To calc nodal mass
  //CalcNodalMassFromVol(); //Repla

  double dt = 1.0e-3; 
  
  Time = 0.0;
  int step_count = 0;
  dt = 1.0;
  double tout = 0;
  
  bool remesh_ = false;
  
  double u_accum[m_dim*m_node_count];
  double u_count[m_dim*m_node_count];
  double relax = 0.5;

  contforce[m_dim*3+2] = -1.0e3;
  
    printf("NODES\n");
    for (int i=0;i<m_node_count;i++){
      for (int d=0;d<m_dim;d++)
        printf("%.4e ", x[m_dim*i+d]);
      printf("\n");
    }  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Main Loop----"<<endl;
  //while (Time < end_t) {
      
  if (step_count % 10000 == 0)
    printf("Step %d, Time %f\n",step_count, Time);  

  
    for (int n=0;n<m_node_count;n++)
      x_old[n]=x[n];  


    // For each element
    cout << "1. Calculating F "<<endl;
    for (int e = 0; e < m_elem_count; e++) {
        // 1) Compute deformation gradient F
        Matrix F(m_dim, m_dim); // 3x3 zero initialized
        Matrix F_old(m_dim, m_dim); // 3x3 zero initialized
        for (int n = 0; n < m_nodxelem; n++) {
            Matrix X(m_dim, 1);
            Matrix X_old(m_dim, 1);
            Matrix gradN(1, m_dim);

            for (int d = 0; d < m_dim; d++) {
                // Current nodal position (x)
                X.Set(d, 0, x[e * m_nodxelem * m_dim + n * m_dim + d]); // careful indexing!
                X_old.Set(d, 0, x_old[e * m_nodxelem * m_dim + n * m_dim + d]); // careful indexing!
                // Shape function gradient in current config
                gradN.Set(0, d, getDerivative(e, 0, d, n)); 
            }
            F += MatMul(X, gradN); // F += x_a ⊗ ∇N_a
            F_old += MatMul(X_old, gradN); // F += x_a ⊗ ∇N_a
        }
        
        cout << "F "<<endl;
        F.Print();
        
        ///// DERIVATION FOR L IS MORE COMMON
        // // 2) Compute Left Cauchy-Green tensor b = F * F^T
        // Matrix b = MatMul(F, F.Transpose());

        // // 3) Compute Almansi strain: e = 0.5 * (I - b^{-1})
        // Matrix I = Identity(m_dim);
        // Matrix b_inv = b.Inv();
        // Matrix e_almansi = (I - b_inv) * 0.5; //almansi
        
        cout << "Calculating Eps "<<endl;
        Matrix eps(m_dim,m_dim) ;
        if (step_count == 0) {
            //eps = Matrix(m_dim, m_dim); // cero
        } else {
            Matrix L = 1.0/dt * (MatMul(F, F_old.Inv()) - Identity(m_dim));
            eps = 0.5 * (L + L.getTranspose());
        }

        cout << "Calculating Strain Voight Notation "<<endl;
        // 4) Convert strain tensor e_almansi (3x3) to 6x1 Voigt vector (engineering strains)
        Matrix strain_voigt(6, 1);
        strain_voigt.Set(0, 0, eps.getVal(0, 0));  // ε_xx
        strain_voigt.Set(1, 0, eps.getVal(1, 1));  // ε_yy
        strain_voigt.Set(2, 0, eps.getVal(2, 2));  // ε_zz
        strain_voigt.Set(3, 0, 2 * eps.getVal(0, 1)); // γ_xy = 2 * ε_xy
        strain_voigt.Set(4, 0, 2 * eps.getVal(1, 2)); // γ_yz
        strain_voigt.Set(5, 0, 2 * eps.getVal(2, 0)); // γ_zx
        
        cout << "2. Calculating D "<<endl;
        // 5) Compute stress σ = D * ε
        Matrix D(6,6);
        D =  mat[e]->getElasticMatrix();
        Matrix stress_voigt = MatMul(D, strain_voigt); // D is 6x1 elastic stiffness matrix
        cout << "Calculating B "<<endl;
        
        Matrix sigma = MatMul(D,eps); /// USE VOIGTH TO MARIX
        //Matrix sigma(3,3);
        //sigma = VoigtToMatrix(stress_voigt);
        
        // 6) Build B matrix (strain-displacement) for the element
        Matrix B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)
        cout <<"Done."<<endl;
        cout << "B mat "<<endl;
        B.Print();
        cout << "m_dim "<<m_dim<<endl;
        // 7) Compute internal force: fint = V_e * B^T * σ
        cout << "Computing internal force"<<endl;
        Matrix fint = MatMul(B.getTranspose(), stress_voigt); //DO NOT TRANSPOSE B DEFITELY
        fint = fint * vol[e];
        cout << "Calculating Kmat "<<endl;
        // // 8.1) Compute tangent stiffness matrix Ktan = V_e * B^T * D * B
        Matrix Kmat = MatMul(B.getTranspose(), MatMul(D, B));
        Kmat = Kmat * vol[e];
        // Kmat = Kmat * vol[e];
        cout << "Kmat "<<endl;
        Kmat.Print();
        
        // // 8.2) (Optional) Compute geometric stiffness Kgeo if needed


        double Ve = vol[e]; // Current volume (updated Lagrangian)

        // 2. Initialize Kgeo (12x12 for 4-node tetrahedron)
        //Matrix& Kgeo = *(m_Kgeo[e]);
        Matrix Kgeo(m_dim*m_node_count,m_dim*m_node_count);
        Kgeo.SetZero();

        // 3. Loop over node pairs (a, b)
        for (int a = 0; a < 4; ++a) {
          // ∇Nᵃ in current config (∂Nᵃ/∂x, ∂Nᵃ/∂y, ∂Nᵃ/∂z)
          Matrix grad_a(3, 1);
          grad_a.Set(0, 0, getDerivative(e, 0, 0, a)); // ∂N/∂x
          grad_a.Set(1, 0, getDerivative(e, 0, 1, a)); // ∂N/∂y
          grad_a.Set(2, 0, getDerivative(e, 0, 2, a)); // ∂N/∂z

          for (int b = 0; b < 4; ++b) {
            // ∇Nᵇ in current config
            Matrix grad_b(3, 1);
            grad_b.Set(0, 0, getDerivative(e, 0, 0, b));
            grad_b.Set(1, 0, getDerivative(e, 0, 1, b));
            grad_b.Set(2, 0, getDerivative(e, 0, 2, b));

            // Compute K_geo(a,b) = (∇Nᵃ)ᵀ · σ · ∇Nᵇ * Ve
            Matrix sigma_grad_b = MatMul(sigma, grad_b); // σ · ∇Nᵇ (3x1)
            Matrix kab = MatMul(grad_a.getTranspose(), sigma_grad_b); // 1x1 scalar
            double k_ab = kab.getVal(0, 0) * Ve;

            // Fill 3x3 block (assumes 3 DOF per node)
            for (int i = 0; i < 3; ++i) {
              Kgeo.Set(3*a + i, 3*b + i, Kgeo.getVal(3*a + i, 3*b + i) + k_ab);
            }
          }
        }
        
      



        // Apply BCs BEFORE the element loop by modifying the element matrices
        for (int e = 0; e < m_elem_count; e++) {
            // Get the element stiffness matrix (Kmat) and force vector (fint)
            //Matrix& K_e = m_Kmat[e];  // Reference to element stiffness
            //Matrix& f_e = fint[e];    // Reference to element force vector

            // Loop over all BC directions (x, y, z)
            for (int dim = 0; dim < m_dim; dim++) {
                // Loop over all BC nodes in this direction
                for (int n = 0; n < bc_count[dim]; n++) {
                    // Get the constrained node and DOF
                    int node = (dim == 0) ? bcx_nod[n] : 
                              ((dim == 1) ? bcy_nod[n] : bcz_nod[n]);
                    int constrained_dof_global = node * m_dim + dim;

                    // Check if this element contains the constrained node
                    for (int ne = 0; ne < m_nodxelem; ne++) {
                        int elem_node = getElemNode(e, ne);
                        if (elem_node == node) {
                            // Local DOF in the element matrix
                            int constrained_dof_local = ne * m_dim + dim;
                            double bc_value = (dim == 0) ? bcx_val[n] : 
                                             ((dim == 1) ? bcy_val[n] : bcz_val[n]);
                                             
                            // --- Modify Kmat ---
                            // Zero out the row and column (except diagonal)
                            for (int k = 0; k < m_nodxelem * m_dim; k++) {
                                if (k == constrained_dof_local) {
                                    Kmat.Set(constrained_dof_local, k, 1.0); // diagonal 1
                                    Kgeo.Set(constrained_dof_local, k, 0.0);
                                } else {
                                    // Primero leemos el valor original antes de modificar
                                    double kij = Kmat.getVal(k, constrained_dof_local);
                                    
                                    // Ajustamos el RHS para mantener el equilibrio
                                    fint.Set(k, 0, fint.getVal(k, 0) - kij * bc_value);

                                    // Ahora anulamos fila y columna
                                    Kmat.Set(constrained_dof_local, k, 0.0); // fila
                                    Kmat.Set(k, constrained_dof_local, 0.0); // columna

                                    Kgeo.Set(constrained_dof_local, k, 0.0);
                                    Kgeo.Set(k, constrained_dof_local, 0.0);
                                }
                            }

                            // Finalmente fijamos el valor impuesto en RHS
                            fint.Set(constrained_dof_local, 0, bc_value);
                            // --- Modify fint ---
                            // Set force to prescribed displacement (if any)
                            

                            cout << "DOF "<< constrained_dof_local<<"BCVAL "<<bc_value<<endl;
                            fint.Set(constrained_dof_local, 0, bc_value);
                        }
                    }
                }
            }
        }
        cout << "Kmat reduced"<<endl;
        Kmat.Print();
        cout << "F"<<endl;
        fint.Print();
        // // 9) Local System: Kmat * Δu_e = fint
        //Matrix delta_u_e = MatMul(Kmat.Inv(),fint) * (-1.0);  // Δu = -K⁻¹·fint (¡signo importante!)
        
        //Matrix du = SolveWithEigen(Kmat, fint);
        Matrix delta_u_e = Matrix::SolveLU( Kmat,fint);
        
        
        cout << "Delta U "<<endl;
        delta_u_e.Print();
        
        // // 10) Distribute Δu_e to nodes
        for (int a = 0; a < m_nodxelem; ++a) {
            for (int d = 0; d < m_dim; ++d) {
                int idx_local = a * m_dim + d;
                //int idx_global = elem_to_node[e][a] * m_dim + d;
                int idx_global = getElemNode(e,a) * m_dim + d;

                u[idx_global] += relax * delta_u_e.getVal(idx_local, 0); // acumulás con relajación
                u_count[idx_global] += 1; // acumulás contribuciones
            }
        }
    
    } // end element loop
  //////////
  
  
  
  ////// BC application.
  for (int dim=0;dim<m_dim;dim++){
    par_loop (n,bc_count[dim]){
      double val;
      //printf("thread %d, Imposing Vel in dim %d, %d Conditions, val %f\n", n, dim, bc_count[dim], bcx_val[n]);
      //printf("BCV dim %d\n", dim);
      // printf("VEL BC \n");

      if (dim == 0)       {
        int ind = m_dim*bcx_nod[n]+dim; 
        u[ind] = bcx_val[n]; 
        u_count[ind]=1;}
      else if (dim == 1)  {
        int ind = m_dim*bcy_nod[n]+dim;
        u[ind] = bcy_val[n];}
      else if (dim == 2)  {
        int ind = m_dim*bcz_nod[n]+dim;
        u[ind] = bcz_val[n]; 
        u_count[ind]=1;}
    }// Node loop
  }//dim loop
  
    // //~ // 11) Average &  actualize nodal pos
    for (int i = 0; i < m_node_count * m_dim; ++i) {
        if (u_count[i] > 0) {
            double delta = u[i] / u_count[i];
            x[i] += delta; // actualizás posición
        }
    }
    
    printf("DISPLACEMENTS\n");
    for (int i=0;i<m_node_count;i++){
      for (int d=0;d<m_dim;d++)
        printf("%.4e ", u[m_dim*i+d]);
      printf("\n");
    }

    // 12) Volver a calcular deformaciones con x actualizado

        // ...





  cout << "Writing output "<<endl;

  // //cout << "Writing output"<<endl;
  // VTKWriter writer2(this, "out.vtk");
  // writer2.writeFile();

  // cout << "Done."<<endl;

  // ///// IF REMESH
  // /////#####################
  // //remesh.WriteDomain();
  // //calcElemJAndDerivatives();    
  
  // //////////////////////////////////////
  
  // VTKWriter writer3(this, "out_remesh.vtk");
  // writer3.writeFile();
  
  // //AFTER WRITE

  // timer.stop();
  // std::cout << "Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  
  }//SOLVE
  
 
//PROBLEMWITH INHERITANCE
 
void host_ Domain_d::ElasticSolve(){
  WallTimer timer;

  AssignMatAddress();
  calcElemJAndDerivatives();
  CalcElemVol();

  //IMPLICIT 
  CalcMaterialStiffElementMatrix();

  par_loop(e, m_elem_count){
    
  }
  
  Solver_Eigen *solver = new Solver_Eigen();
  m_solver = solver;
  m_solver->setDomain(this);
  m_solver->Allocate();
  cout << "Assemblying matrix "<<endl;
  m_solver->assemblyGlobalMatrix();

  
  m_solver->applyDirichletBCs();
  cout << "Solving system"<<endl;
  m_solver->Solve();
  
  
  }//ELASTICSOLVE
  

    
};
