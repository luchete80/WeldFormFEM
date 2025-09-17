#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"
#include "Matrix.h"

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>  // ofstream
#include <iomanip>

#ifdef BUILD_REMESH
#include "ReMesher.h"
#endif

#include <omp.h>

#include "Solver_Eigen.h"
#include "tensor3.C"

using namespace std;




namespace MetFEM{
  
//~ void Domain_d::ImposeBCU(int dim) {
    //~ int* bc_nodes = nullptr;
    //~ double* bc_values = nullptr;
    
    //~ // Seleccionar el array correcto según la dimensión
    //~ switch (dim) {
        //~ case 0: bc_nodes = bcx_nod; bc_values = bcx_val; break;
        //~ case 1: bc_nodes = bcy_nod; bc_values = bcy_val; break;
        //~ case 2: bc_nodes = bcz_nod; bc_values = bcz_val; break;
        //~ default: return;
    //~ }
    
    //~ // Aplicar BCs de desplazamiento de forma INCREMENTAL
    //~ //par_loop (i, bc_count[dim]) {
      //~ for (int i=0;i<bc_count[dim];i++){
        //~ int node = bc_nodes[i];
        //~ double target_disp = bc_values[i];
        //~ double current_disp = u[node * m_dim + dim];
        
        //~ // Calcular el incremento necesario
        //~ double delta_disp = target_disp - current_disp;
        
        //~ // Aplicar como condición de contorno
        //~ int gdof = node * m_dim + dim;
        //~ cout << "gdof "<<gdof <<", dim "<<dim <<endl;
        //~ m_solver->setDirichletBC(gdof, delta_disp);
    //~ }
//~ }

    void Domain_d::ImposeBCU(int dim/*, double load_factor*/) {
        int* bc_nodes = nullptr;
        double* current_bc_values = nullptr;
        std::vector<double>* original_values = nullptr;
        
        // Seleccionar arrays según dimensión
        switch (dim) {
            case 0: 
                bc_nodes = bcx_nod; 
                current_bc_values = bcx_val;
                original_values = &original_bcx_val;
                break;
            case 1: 
                bc_nodes = bcy_nod; 
                current_bc_values = bcy_val;
                original_values = &original_bcy_val;
                break;
            case 2: 
                bc_nodes = bcz_nod; 
                current_bc_values = bcz_val;
                original_values = &original_bcz_val;
                break;
            default: return;
        }
        
        // Aplicar BCs de forma incremental según load_factor
        par_loop (i, bc_count[dim]) {
            double target_disp = (*original_values)[i] /** load_factor*/;
            double current_disp = u[bc_nodes[i] * m_dim + dim];
            double delta_disp = target_disp - current_disp;
            
            //cout << "Writing bc "<< i <<", node "<<bc_nodes[i]<<endl;
            // Sobrescribir el valor en el array de BCs
            current_bc_values[i] = delta_disp; // ← Ahora guarda el incremento

        }
    }



void Domain_d::SolveStaticDisplacement() {  
  WallTimer timer;
  std::ofstream of("Contact_Forces.csv", std::ios::out);
  
  // 1. ELIMINAR VARIABLES DINÁMICAS
  // Remover: prev_v, prev_a, v, a, m_mdiag (masas)
  // Mantener: x, u (desplazamientos), f (fuerzas)

  double *x_initial;
  #ifndef BUILD_GPU
    x_initial = new double [m_dim * m_node_count];
  #endif
  memcpy(x_initial, x, sizeof(double) * m_node_count * m_dim);

  double* delta_u = new double[m_dim * m_node_count];
  memset(delta_u, 0, sizeof(double) * m_dim * m_node_count);
          
  //~ cout << "Initializing Values ..." << endl;
  //~ //InitValuesStatic(); // Nueva función inicialización estática

  //~ cout << "Done." << endl;
  setOriginalBCs(); //For incremental
    
  cout << "Imposing bcs"<<endl;
  // 2. ELIMINAR INICIALIZACIÓN DE VELOCIDADES/ACELERACIONES
  for (int d = 0; d < m_dim; d++) {
    ImposeBCU(d); // Nueva función para BCs de desplazamiento
  }
  cout << "Done"<<endl;

  // 3. CONFIGURACIÓN SOLVER ESTÁTICO
  Solver_Eigen* solver = new Solver_Eigen();
  m_solver = solver;
  m_solver->setDomain(this);
  cout << "Allocating "<<endl;
  m_solver->Allocate();
  cout << "Done."<<endl;
  

  double* u_previous = new double[m_dim * m_node_count]; // Desplazamiento anterior
  
  
  // 4. ELIMINAR TODO LO RELACIONADO CON TIEMPO
  // Remover: Time, dt, step_count, etc.
  int load_step = 0;
  const int max_load_steps = 1;
  double load_factor = 0.0;
  const double max_load_factor = 1.0;
  
    ////IMPLICIT DEFS (SAVINGS MEM)
  #ifndef BUILD_GPU
    std::vector<Matrix> Bmat_per_thread(Nproc);
    std::vector<Matrix> sig_per_thread(Nproc);
  #else
    
  #endif
  
  
  // 5. LOOP PRINCIPAL DE CARGA (no de tiempo)
  while (load_factor < max_load_factor && load_step < max_load_steps) {
    
    load_factor += 1.0 / max_load_steps;
    //ApplyLoads(load_factor); // Aplicar cargas incrementalmente

    // 6. NEWTON-RAPHSON PARA NO LINEALIDAD GEOMÉTRICA/MATERIAL
    bool converged = false;
    int max_iter = 20;
    double tolerance = 1e-6;
    
    for (int iter = 0; iter < max_iter && !converged; iter++) {
      
      // 7. ACTUALIZAR GEOMETRÍA
      for (int i = 0; i < m_node_count * m_dim; i++) {
        x[i] = x_initial[i] + u[i]; // Posición actualizada
      }

      // 8. CALCULAR DEFORMACIONES Y TENSIONES
      calcElemJAndDerivatives();
      CalcElemVol();
      
      if (m_press_algorithm == 0)
        calcElemPressure();
      else if (m_press_algorithm == 1)
        calcElemPressureANP();

      CalcStressStrain(0.0); // dt = 0 para estático

      // 9. CALCULAR FUERZAS INTERNAS
      calcElemForces();
      calcElemHourglassForces();
      
      if (contact)
        CalcContactForces();


      ///// CLASSSIC (OLD) ELASTIC APPROACH
      CalcMaterialStiffElementMatrix();
      //~ solver->assemblyGlobalMatrix();

    
      // 10. ENSAMBLAR SISTEMA
      solver->setZero();
      solver->beginAssembly();

    
    
    cout << "Building matrices "<<endl;
    /////////////////////// THIS IS BEB
    //par_loop(e,m_elem_count){
    for (int e=0;e<m_elem_count;e++){
          //cout << "Element "<<e<<endl;
          
          // 6) Build B matrix (strain-displacement) for the element
          //int tid = omp_get_thread_num();

          //Matrix &B = Bmat_per_thread[tid];
          Matrix B;
          //// HERE B is in fact BxdetJ
          B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)
          B = B *(1.0/m_detJ[e]);

          //cout << "B mat "<<endl;
          //B.Print();

          // 7) Compute internal force: fint = V_e * B^T * σ

          Matrix stress_voigt = FlatSymToVoigt(m_sigma,m_dim,m_nodxelem);
          //CHANGE TO FORCES TO MATRIX! already calculated
          Matrix fint = MatMul(B.getTranspose(), stress_voigt); //DO NOT TRANSPOSE B DEFITELY
 
          fint = fint * vol[e];

          ///TANGENT!
          // // 8.1) Compute tangent stiffness matrix Ktan = V_e * B^T * D * B
          Matrix D(6,6);
          D =  mat[e]->getElasticMatrix();
          
          Matrix Kmat = MatMul(B.getTranspose(), MatMul(D, B));
          Kmat = Kmat * (1.0/6.0*m_detJ[e]); // B is B x detJ

          double Ve = vol[e]; // Current volume (updated Lagrangian)

          ///////////////////////////////////////////////////
          /////////// IMPORTANT!!! --A LOT-- FASTER (LESS PRODUCTS) THAN: Kgeo = G^T sigma G
          // 2. Initialize Kgeo (12x12 for 4-node tetrahedron)
          //Matrix& Kgeo = *(m_Kgeo[e]);
          //~ Matrix Kgeo(m_dim*m_nodxelem,m_dim*m_nodxelem);
          //~ Kgeo.SetZero();
          
          //~ // // 3. Loop over node pairs (a, b)
          //~ // // REMEMBER DERIVATIVES ARE AFFECTED BY DETJ
          //~ for (int a = 0; a < 4; ++a) {
            //~ // ∇Nᵃ in current config (∂Nᵃ/∂x, ∂Nᵃ/∂y, ∂Nᵃ/∂z)
            //~ Matrix grad_a(3, 1);
            //~ grad_a.Set(0, 0, getDerivative(e, 0, 0, a)); // ∂N/∂x
            //~ grad_a.Set(1, 0, getDerivative(e, 0, 1, a)); // ∂N/∂y
            //~ grad_a.Set(2, 0, getDerivative(e, 0, 2, a)); // ∂N/∂z

            //~ for (int b = 0; b < 4; ++b) {
              //~ // ∇Nᵇ in current config
              //~ Matrix grad_b(3, 1);
              //~ grad_b.Set(0, 0, getDerivative(e, 0, 0, b));
              //~ grad_b.Set(1, 0, getDerivative(e, 0, 1, b));
              //~ grad_b.Set(2, 0, getDerivative(e, 0, 2, b));

              //~ // Compute K_geo(a,b) = (∇Nᵃ)ᵀ · σ · ∇Nᵇ * Ve
              //~ Matrix sigma_grad_b = MatMul(FlatSymToMatrix(m_sigma), grad_b); // σ · ∇Nᵇ (3x1)
              //~ Matrix kab = MatMul(grad_a.getTranspose(), sigma_grad_b); // 1x1 scalar
              //~ double k_ab = kab.getVal(0, 0) * Ve;

              //~ // Fill 3x3 block (assumes 3 DOF per node)
              //~ for (int i = 0; i < 3; ++i) {
                //~ Kgeo.Set(3*a + i, 3*b + i, Kgeo.getVal(3*a + i, 3*b + i) + k_ab);
              //~ }
            //~ }
          //~ }

          //~ Kgeo = Kgeo * (1.0/(6.0*m_detJ[e]));
          //~ Matrix K = Kgeo + Kmat;
          
          Matrix K = Kmat;
          //cout << "K MATRIX"<<endl;
          //K.Print();
          //cout << "Ke "<<endl;
          //m_Kmat[e]->Print();
          

          Matrix R(m_dim*m_nodxelem,1);
          for (int i = 0; i < m_nodxelem; i++) {
            //int node = getElemNode(e, i % m_nodxelem);
            for (int d=0;d<m_dim;d++){
            //cout << "NODE, DIM "<<i<<","<<d<<", fint mat"<<fint.getVal(m_dim*i+d,0)<<", fel "<<m_f_elem[i*m_dim+d]<<endl;
            R.Set(i,0,0.0); //ADD EXTERNAL ELEMENT FORCES
            //R.Set(m_dim*i+d,0,-m_f_elem[i*m_dim+d]); //ADD EXTERNAL ELEMENT FORCES
            }
          }
          solver->assembleElement(e, K);
          //solver->assembleResidual(e,R);//SHOULD BE NEGATIVE!  
          
      }//elem
      solver->finalizeAssembly();
      
      
  
      cout << "Done "<<endl;
      
      // Aplicar condiciones de contorno
      for (int d = 0; d < m_dim; d++)
        ImposeBCU(d);
        
      for (int n = 0; n < m_node_count*m_dim; n++)      
        solver->addToR(n,contforce[n]); //EXTERNAL FORCES
        
      m_solver->applyDirichletBCs();
      cout << "Solving "<<endl;
      // 11. RESOLVER SISTEMA
      m_solver->Solve();
      
      // 12. ACTUALIZAR DESPLAZAMIENTOS
      double max_residual = 0.0;
      for (int n = 0; n < m_node_count; n++) {
        for (int d = 0; d < m_dim; d++) {
          int idx = n * m_dim + d;
          double du = m_solver->getU(n, d);
          delta_u[n * m_dim + d] += du;
          //u[idx] += du;
          max_residual = std::max(max_residual, std::abs(du));
        }
      }
      
      // 13. VERIFICAR CONVERGENCIA
      if (max_residual < tolerance) {
        converged = true;
        cout << "Converged in " << iter + 1 << " iterations" << endl;

            for (int i = 0; i < m_node_count * m_dim; i++) {
                u[i] += delta_u[i];
            }
            cout << "Converged in " << iter + 1 << " iterations" << endl;
            
            // 14. CALCULAR TENSIONES FINALES
            for (int i = 0; i < m_node_count * m_dim; i++) {
                x[i] = x_initial[i] + u[i];
            }
            
            calcElemJAndDerivatives();
            CalcElemVol();
            CalcStressStrain(0.0);
                  
      
      }/////NR ITER 
    
    
    }////// NR ITER
    
    // 14. OUTPUT Y REMALLADO
    if (load_step % 10 == 0) {
      string outfname = "out_load_" + std::to_string(load_factor) + ".vtk";
      VTKWriter writer(this, outfname.c_str());
      writer.writeFile();
    }
    
    load_step++;
  }

  // 15. LIMPIEZA
  delete[] delta_u;
  delete[] u_previous;
  of.close();
}


    
};


