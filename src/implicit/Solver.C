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
  
//~ void radialReturnJ2(const Matrix& D, double dt, 
                 //~ const Matrix& sigma_old, double ep_old,
                 //~ const MaterialProps& props,
                 //~ Matrix& sigma_new, double& ep_new,
                 //~ Matrix& C_ep) {
  
  //~ // Elastic predictor
  //~ Matrix deps = D * dt;
  //~ Matrix sigma_trial = sigma_old + 2*props.G*deviatoric(deps) 
                     //~ + props.K*trace(deps)*Matrix::Identity();
  
  //~ // Check yield condition
  //~ double f_trial = equivalentStress(sigma_trial) 
                 //~ - props.yieldStress(ep_old);
  
  //~ if (f_trial <= 0) {
      //~ // Elastic step
      //~ sigma_new = sigma_trial;
      //~ ep_new = ep_old;
      //~ C_ep = elasticTangent(props);
  //~ } else {
      //~ // Plastic correction
      //~ double delta_gamma = f_trial / (3*props.G + props.H);
      //~ Matrix n = deviatoric(sigma_trial) / equivalentStress(sigma_trial);
      //~ sigma_new = sigma_trial - 2*props.G*delta_gamma*n;
      //~ ep_new = ep_old + sqrt(2.0/3.0)*delta_gamma;
      
      //~ // Consistent tangent
      //~ C_ep = consistentTangent(props, sigma_trial, delta_gamma);
  //~ }
//~ }

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
  dt = 1.0;
  double tout = 0;
  
  bool remesh_ = false;
  
  double u_accum[m_dim*m_node_count];
  double u_count[m_dim*m_node_count];
  double relax = 0.5;

  
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
        Matrix stress_voigt = MatMul(D, strain_voigt); // D is 6x6 elastic stiffness matrix
        cout << "Calculating B "<<endl;
        
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
        // Kmat = Kmat * vol[e];
        cout << "Kmat "<<endl;
        Kmat.Print();
        
        // // 8.2) (Optional) Compute geometric stiffness Kgeo if needed

        // // 9) Local System: Kmat * Δu_e = fint
        Matrix delta_u_e = MatMul(Kmat.Inv(),fint) * (-1.0);  // Δu = -K⁻¹·fint (¡signo importante!)
        
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

  for (int dim=0;dim<m_dim;dim++){
    par_loop (n,bc_count[dim]){
      double val;
      //printf("thread %d, Imposing Vel in dim %d, %d Conditions, val %f\n", n, dim, bc_count[dim], bcx_val[n]);
      //printf("BCV dim %d\n", dim);
      // printf("VEL BC \n");
      if (dim == 0)       {/*printf ("dim %d node %f, val %d\n",dim,bcx_nod[n],bcx_val[n]); */ u[m_dim*bcx_nod[n]+dim] = bcx_val[n]; }
      else if (dim == 1)  {/*printf ("dim %d node %d val %f \n",dim,bcy_nod[n], bcy_val[n]);*/ u[m_dim*bcy_nod[n]+dim] = bcy_val[n];}
      else if (dim == 2)  {/*printf ("dim %d node %f, val %d\n",dim,bcz_nod[n],bcz_val[n]);*/  u[m_dim*bcz_nod[n]+dim] = bcz_val[n]; }
    }
  }
    // //~ // 11) Average &  actualize nodal pos
    for (int i = 0; i < m_node_count * m_dim; ++i) {
        if (u_count[i] > 0) {
            double delta = u[i] / u_count[i];
            x[i] += delta; // actualizás posición
        }
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
