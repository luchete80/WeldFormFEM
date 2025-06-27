#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"

#include "Matrix.h"

#include "Solver_Eigen.h"

// #ifdef BUILD_REMESH
// #include "ReMesher.h"
// #endif

using namespace std;

namespace MetFEM{

void host_ Domain_d::Solve(){
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
  double tout = 0;
  
  bool remesh_ = false;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Main Loop----"<<endl;
  //while (Time < end_t) {
      
  if (step_count % 10000 == 0)
    printf("Step %d, Time %f\n",step_count, Time);  

  
  


// For each element
for (int e = 0; e < m_elem_count; e++) {
    // 1) Compute deformation gradient F
    Matrix F(m_dim, m_dim); // 3x3 zero initialized
    for (int n = 0; n < m_nodxelem; n++) {
        Matrix X(m_dim, 1);
        Matrix gradN(1, m_dim);

        for (int d = 0; d < m_dim; d++) {
            // Current nodal position (x)
            X.Set(d, 0, x[e * m_nodxelem * m_dim + n * m_dim + d]); // careful indexing!
            // Shape function gradient in current config
            gradN.Set(0, d, getDerivative(e, 0, d, n)); 
        }
        F += MatMul(X, gradN); // F += x_a ⊗ ∇N_a
    }

    // 2) Compute Left Cauchy-Green tensor b = F * F^T
    Matrix b = MatMul(F, F.Transpose());

    // 3) Compute Almansi strain: e = 0.5 * (I - b^{-1})
    Matrix I = Identity(m_dim);
    Matrix b_inv = b.Inv();
    Matrix e_almansi = (I - b_inv) * 0.5;

    // 4) Convert strain tensor e_almansi (3x3) to 6x1 Voigt vector (engineering strains)
    Matrix strain_voigt(6, 1);
    strain_voigt.Set(0, 0, e_almansi.getVal(0, 0));  // ε_xx
    strain_voigt.Set(1, 0, e_almansi.getVal(1, 1));  // ε_yy
    strain_voigt.Set(2, 0, e_almansi.getVal(2, 2));  // ε_zz
    strain_voigt.Set(3, 0, 2 * e_almansi.getVal(0, 1)); // γ_xy = 2 * ε_xy
    strain_voigt.Set(4, 0, 2 * e_almansi.getVal(1, 2)); // γ_yz
    strain_voigt.Set(5, 0, 2 * e_almansi.getVal(2, 0)); // γ_zx

    // 5) Compute stress σ = D * ε
    Matrix D(6,6);
    D =  mat[e]->getElasticMatrix();
    Matrix stress_voigt = MatMul(D, strain_voigt); // D is 6x6 elastic stiffness matrix

    // 6) Build B matrix (strain-displacement) for the element
    Matrix B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)

    // 7) Compute internal force: fint = V_e * B^T * σ
    Matrix fint = MatMul(B.Transpose(), stress_voigt);
    fint = fint * vol[e];

    // 8) Compute tangent stiffness matrix Ktan = V_e * B^T * D * B
    Matrix Kmat = MatMul(B.Transpose(), MatMul(D, B));
    Kmat = Kmat * vol[e];

    // 9) (Optional) Compute geometric stiffness Kgeo if needed

    // 10) Assemble fint and Kmat into global system (your assembler)

    // ...

} // end element loop



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
