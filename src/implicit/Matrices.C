/*************************************************************************/
/*  Matrices.C                                                   */
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
#include "Matrix.h"

#include <iostream>
#include <unordered_map>  // For triplet storage

#include "Tensor.h"

using namespace std;

namespace MetFEM{
  
  
//CALCULATE DEFORMATION GRADIENT F
/// F = I + Grad * dU


/// ALMANSI
//// e = 1/2(I-FT F-1)

 //////[dN1/dx    0      0 ]
     //[   0   dN1/dy    0]
     //[   0      0    dN1/dx ]    
     //[   dNdy  dN/dx    0 ]    

// Matrix Domain_d::getB(int &e){
  
      // for (int i=0;i<m_nodxelem;i++)
      // for (int d=0;d<m_dim;d++)
        // Bmat.Set(d, m_dim*i+d, getDerivative(e, 0, d, m_dim*i));
      
    // if (m_dim==3){
      // for (int i=0;i<m_nodxelem;i++)
        // for (int d=0;d<m_dim;d++){
          // int k = d+1;if (k==m_dim) k = 0;
            // // d/dy d/dx 0   
            // printf("i %d j %d der %d\n",m_dim+d,m_dim*i+d, k);
            // printf("i %d j %d der %d\n",m_dim+d,m_dim*i+k, d);
            // Bmat.Set(m_dim+d, m_dim*i+d, getDerivative(e, 0, k, i));
            // Bmat.Set(m_dim+d, m_dim*i+k,getDerivative(e, 0, d, i));
        // }
    // }
// }

///// USES VOL OR DETERMINANT
/////// IN CUDA ISPREFERRED NOT TO SAVE IT
Matrix Domain_d::getElemBMatrix(const int &e){
  Matrix B(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?
    if (m_dim == 3){
    //3D Voigt notation, where:
    for (int i=0;i<m_nodxelem;i++){
      int base = m_dim * i;
      double dN_dx = getDerivative(e, 0, 0, i);
      double dN_dy = getDerivative(e, 0, 1, i);
      double dN_dz = getDerivative(e, 0, 2, i);
      
      /// BATHE'S ORDER, FOR tau exx eyy ezz exy eyz ezx
      //cout << " dN_dx dN_dy dN_dz: %f %f %f "<<dN_dx<<", "<<dN_dy<<", "<<dN_dz<<endl;
      B.Set(0, base + 0, dN_dx);  // ε_xx
      B.Set(1, base + 1, dN_dy);  // ε_yy
      B.Set(2, base + 2, dN_dz);  // ε_zz
      
      B.Set(3, base + 0, dN_dy);  // ε_xy
      B.Set(3, base + 1, dN_dx);
      B.Set(4, base + 1, dN_dz);  // ε_yz
      B.Set(4, base + 2, dN_dy);
      B.Set(5, base + 2, dN_dx);  // ε_zx
      B.Set(5, base + 0, dN_dz);
    }
  }
  return B;
}


void Domain_d::CalcMaterialStiffElementMatrix(){
  par_loop(e, m_elem_count){
    /// TODO: CHANGE FOR DIM = 2
    Matrix B(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?
    if (m_dim == 3){
    //3D Voigt notation, where:
    for (int i=0;i<m_nodxelem;i++){
      int base = m_dim * i;
      double dN_dx = getDerivative(e, 0, 0, i);
      double dN_dy = getDerivative(e, 0, 1, i);
      double dN_dz = getDerivative(e, 0, 2, i);
      
      //cout << " dN_dx dN_dy dN_dz: %f %f %f "<<dN_dx<<", "<<dN_dy<<", "<<dN_dz<<endl;
      B.Set(0, base + 0, dN_dx);  // ε_xx
      B.Set(1, base + 1, dN_dy);  // ε_yy
      B.Set(2, base + 2, dN_dz);  // ε_zz
      
      B.Set(3, base + 0, dN_dy);  // ε_xy
      B.Set(3, base + 1, dN_dx);
      B.Set(4, base + 1, dN_dz);  // ε_yz
      B.Set(4, base + 2, dN_dy);
      B.Set(5, base + 2, dN_dx);  // ε_zx
      B.Set(5, base + 0, dN_dz);
    }
  }
    B = B *(1.0/m_detJ[e]);
    
    //printf ("BMAT\n");
    //B.Print();
    Matrix BT(m_nodxelem* m_dim,2*m_dim);
    BT=B.getTranspose();
    Matrix D(6,6);
    
    double G  = mat[e]->Elastic().G();
    //cout << "Material G "<<G<<endl;
    double E  = mat[e]->Elastic().E();
    double nu = mat[e]->Elastic().Poisson();
    double lambda  = (E * nu) /((1.0+nu)*(1.0-2.0*nu)); 
    D.Set(0,1, lambda);               D.Set(0,2, lambda);
    D.Set(1,0, lambda);               D.Set(1,2, lambda);
    D.Set(2,0, lambda);               D.Set(2,1, lambda);
    
    //printf("Cmat\n");
    
    for (int d=0;d<3;d++) D.Set(d,d,lambda+2.0*G);
    for (int d=3;d<6;d++) D.Set(d,d,G);    
    //D.Print();
    
    MatMul(MatMul(BT,D),B, m_Kmat[e]);
    //printf("K ELEM\n");
    //m_Kmat[e]->Print();

    // Compute element stiffness: K_e = ∫Bᵀ·D·B dV ≈ Bᵀ·D·B * Ve
    *(m_Kmat[e]) = *(m_Kmat[e]) * vol[e];  // Multiply by element volume (Ve)
    
    //cout << "KMAT"<<endl;
    //(m_Kmat[e])->Print();
    
  }//element
  //cout << "Done."<<endl;
  
}

dev_t Matrix Domain_d::CalcElementStressAndTangent(int e, double dt) {
    Matrix K_elem(m_nodxelem * m_dim, m_nodxelem * m_dim);
    K_elem.SetZero();
    tensor3 Sigma_final, ShearStress_final;

    int gp = 0;
      int offset_s = e * m_gp_count + gp;    // SCALAR offset
      int offset_t = offset_s * 6;           // TENSOR offset
    double G = mat[e]->Elastic().G();
    // 1. PREDICTOR ELÁSTICO (CORRECTO)
    tensor3 ShearStress_old = FromFlatSym(m_tau, offset_t);
    tensor3 StrRate = FromFlatSym(m_str_rate, offset_t);
    tensor3 RotRate = FromFlatAntiSym(m_rot_rate, offset_t);
    
    // Término de Jaumann (CORRECTO)
    tensor3 SRT = ShearStress_old * Trans(RotRate);
    tensor3 RS = RotRate * ShearStress_old;
    tensor3 StrRateDev = StrRate - (1.0/3.0)*Trace(StrRate)*Identity();
    
    tensor3 ShearStress_trial = ShearStress_old + dt * (2.0 * G * StrRateDev + SRT + RS);
    tensor3 Sigma_trial = -p[offset_s] * Identity() + ShearStress_trial;
    
    // 2. CALCULAR DESVIADOR TRIAL (CORRECTO)
    tensor3 s_trial = Sigma_trial - (1.0/3.0)*Trace(Sigma_trial)*Identity();
    //double J2_trial = 0.5 * DoubleContraction(s_trial, s_trial);
    double J2_trial = 0.5 * (s_trial.xx*s_trial.xx + 2.0*s_trial.xy*s_trial.xy + 
                        2.0*s_trial.xz*s_trial.xz + s_trial.yy*s_trial.yy + 
                        2.0*s_trial.yz*s_trial.yz + s_trial.zz*s_trial.zz);
    double sig_trial = sqrt(3.0 * J2_trial);
    
    // 3. TENSION DE FLUENCIA (usar pl_strain_old)
    double sigma_y ;//= CalcHollomonYieldStress(pl_strain_old, mat[e]);
      if (mat[e]->Material_model == HOLLOMON) {
          sigma_y = CalcHollomonYieldStress(pl_strain[e], mat[e]);
      } else if (mat[e]->Material_model == JOHNSON_COOK) {
          //double eff_strain_rate = sqrt(0.5*(  ));
          //sigma_y = CalcJohnsonCookYieldStress(pl_strain[e], eff_strain_rate, T[e], mat[e]);
      } else {
        
        sigma_y = mat[e]->sy0;
        }
    
    Matrix D_gp(6,6);
    // 4. VERIFICAR PLASTICIDAD
    if (sig_trial <= sigma_y) {
        // PASO ELÁSTICO
        Sigma_final = Sigma_trial;
        ShearStress_final = ShearStress_trial;
        D_gp = mat[e]->getElasticMatrix();
    } else {///RADIAL RETURN
        
        // 4.1. Hardening modulus H (TANGENTE, not secant)
        double H = CalcHollomonTangentModulus(pl_strain[e], mat[e]);
        
        // 4.2. Plastic  (FORMULACIÓN CLÁSICA)
        double delta_gamma = (sig_trial - sigma_y) / (3.0 * G + H);
        
        // 4.3. Dirección del flujo plástico
        tensor3 n = s_trial * (1.0 / sig_trial);
        
        // 4.4. RETORNO RADIAL (FORMULACIÓN CLÁSICA)
        tensor3 s_final = s_trial - 2.0 * G * delta_gamma * n;
        
        // 4.5. Actualizar deformación plástica
        double delta_ep = sqrt(2.0/3.0) * delta_gamma;
        double pl_strain_new = pl_strain[e] + delta_ep;
        
        // 4.6. Tensión final
        ShearStress_final = s_final;
        Sigma_final = -p[offset_s] * Identity() + s_final;
        
        // 4.7. Matriz tangente consistente (CORRECTA)
        D_gp = getConsistentPlasticTangentMatrix(s_trial, sig_trial, G, H);
        
        // Guardar para actualización posterior
        //pl_strain_temp[e] = pl_strain_new;
    }

      // // 6. CALCULATE ELEMENT STIFFNESS MATRIX+
      //cout << "calculating"<<endl;
      Matrix B = getElemBMatrix(e);
      B = B *(1.0/m_detJ[e]);

      //Matrix K_gp = MatMul(MatMul(B.Transpose(), D_gp), B) * vol[e] * m_detJ[e];
      Matrix K_gp = MatMul(B.getTranspose(), MatMul(D_gp, B));
      K_elem = K_elem + K_gp;
      K_elem = K_elem * (1.0/6.0*m_detJ[e]); // B is B x detJ     

      // // 7. UPDATE STRESSES (temporary, will be committed if converged)
      ToFlatSymPtr(Sigma_final, m_sigma, offset_t);
      ToFlatSymPtr(ShearStress_final, m_tau, offset_t);
            

    return K_elem;
}


Matrix Domain_d::getConsistentPlasticTangentMatrix(const tensor3 &s_trial,
                                         double sig_trial,
                                         double G,
                                         double H)
{
    // De debe estar en NOTACIÓN INGENIERIL (coincide con tu B)
    Matrix De = mat[0]->getElasticMatrix(); // 6x6

    //------------------------------------------------------
    // 1) n_tensor = (3/2) * s / sigma_eq   (derivada de von Mises)
    //------------------------------------------------------
    tensor3 n_t;
    if (sig_trial > 0.0) {
        double factor = 1.5 / sig_trial;   // = 3/2 * 1/sigma_eq
        n_t = s_trial * factor;
    } else {
        n_t = tensor3(); // cero
    }

    //------------------------------------------------------
    // 2) Convertir a Voigt-ingenieril:
    //    [xx, yy, zz, gamma_xy, gamma_yz, gamma_zx]
    //    gamma = 2 * epsilon
    //------------------------------------------------------
    Matrix n_voigt(6,1);
    n_voigt.Set(0,0, n_t.xx);
    n_voigt.Set(1,0, n_t.yy);
    n_voigt.Set(2,0, n_t.zz);
    
    //Being B defined enginnering [not with 1/2]
    n_voigt.Set(3,0, 2.0 * n_t.xy);
    n_voigt.Set(4,0, 2.0 * n_t.yz);
    n_voigt.Set(5,0, 2.0 * n_t.zx);

    Matrix nT = n_voigt.getTranspose();

    //------------------------------------------------------
    // 3) De * n   → vector 6×1
    //------------------------------------------------------
    Matrix De_n = MatMul(De, n_voigt);

    //------------------------------------------------------
    // 4) denom = H + n^T De n
    //------------------------------------------------------
    Matrix temp = MatMul(nT, De_n);
    double denom = H + temp.getVal(0,0);

    const double tol = 1e-14;
    if (std::abs(denom) < tol)
        denom = (denom >= 0 ? tol : -tol);

    //------------------------------------------------------
    // 5) correction = (De*n)(De*n)^T / denom
    //------------------------------------------------------
    Matrix De_nT = De_n.getTranspose();
    Matrix corr = MatMul(De_n, De_nT) * (1.0 / denom);

    //------------------------------------------------------
    // 6) D_ep = De - correction
    //------------------------------------------------------
    Matrix Dep = De - corr;

    // Simetrizar (para evitar drift numérico)
    for (int i=0; i<6; i++){
        for (int j=i+1; j<6; j++){
            double s = 0.5*(Dep.getVal(i,j)+Dep.getVal(j,i));
            Dep.Set(i,j,s);
            Dep.Set(j,i,s);
        }
    }

    return Dep;
}


/////// MATRIC CALC MAIN FUNCTION
/////// KTG AND KGEO IN ONE PLACE

// #ifndef BUILD_GPU
// // THREAD STORAGE (AVOIDS CPU OVERHEAD!!!) (ej. con OpenMP)
// std::vector<Matrix> Kmat_per_thread(num_threads);
// std::vector<Matrix> Kgeo_per_thread(num_threads);

// #pragma omp parallel for
// for (int e = 0; e < m_elem_count; e++) {
    // int tid = omp_get_thread_num();
    // Matrix& Kmat = Kmat_per_thread[tid]; // Reutiliza memoria
    // Matrix& Kgeo = Kgeo_per_thread[tid];
    
    // computeElementMatrices(e, Kmat, Kgeo); // Rellena matrices
    // // ... resto del cálculo ...
// }
//#else
// Si usas GPU:
// cpp
// // Kernel CUDA/HIP (sin almacenamiento global)
// __global__ void computeElements(Matrix* global_fint, double* u) {
    // int e = blockIdx.x * blockDim.x + threadIdx.x;
    // if (e >= m_elem_count) return;

    // // Matrices locales en registros/memoria compartida
    // Matrix Kmat(24,24); 
    // Matrix Kgeo(24,24);
    
    // // Cálculos
    // computeElementMatrices(e, Kmat, Kgeo, u);
    // // ... scatter a global_fint ...
// }




/*
void Domain_d::CalcGeomStiffElementMatrix() {
  par_loop(e, m_elem_count) {
   
    Matrix& stress = *(m_stress_voigt[e]);  // 3x3 Cauchy stress tensor (assumed full matrix, not Voigt)
    double Ve = vol[e];                     // element volume

    Matrix& Kgeo = *(m_Kgeo[e]);            // 12x12 matrix (4 nodes × 3 DOF)
    Kgeo.SetZero();

    for (int a = 0; a < 4; ++a) {
      // Get ∇Nᵃ = (dN_dx, dN_dy, dN_dz)
      Matrix grad_a(3, 1);
      grad_a.Set(0, 0, getDerivative(e, 0, 0, a)); // ∂N/∂x
      grad_a.Set(1, 0, getDerivative(e, 0, 1, a)); // ∂N/∂y
      grad_a.Set(2, 0, getDerivative(e, 0, 2, a)); // ∂N/∂z

      for (int b = 0; b < 4; ++b) {
        Matrix grad_b(3, 1);
        grad_b.Set(0, 0, getDerivative(e, 0, 0, b));
        grad_b.Set(1, 0, getDerivative(e, 0, 1, b));
        grad_b.Set(2, 0, getDerivative(e, 0, 2, b));

        // Multiply: sigma @ grad_b (3x1)
        Matrix sigma_grad_b = MatMul(stress, grad_b);  // 3x1

        // Then compute: kab = sigma_grad_b @ grad_aᵀ (3x3)
        Matrix grad_a_T = grad_a.getTranspose();       // 1x3
        Matrix kab = MatMul(sigma_grad_b, grad_a_T);   // 3x3

        // Add to Kgeo at block (3a:3a+3, 3b:3b+3)
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            double old_val = Kgeo.getVal(3*a + i, 3*b + j);
            double contrib = Ve * kab.getVal(i, j);
            Kgeo.Set(3*a + i, 3*b + j, old_val + contrib);
          }
        }
      }
    }

    //Kgeo.Print();
  }
}
*/


/// KGEO - NO B APPROACH
//dev_t void Domain_d::CalcGeomStiffElementMatrix(){
   //# Insert Kab into 12x12 Kgeo_local at block [a][b]
 // Matrix Bgeo(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?  
  // # - Kgeo_local: 12x12 zero matrix (4 nodes, 3 DOFs each)

  // for a in range(4):  # loop over nodes
      // for b in range(4):
          // grad_Na = grad_N[a]  # (3x1)
          // grad_Nb = grad_N[b]  # (3x1)
          
          // # Compute 3x3 block: K_ab = Ve * (grad_Na^T * sigma * grad_Nb) * I3
          // # ACTUALLY: K_ab = Ve * outer(grad_Na, grad_Nb) : sigma (contracted)
          // # But your scalar coeff * I3 is a common APPROXIMATION (see notes below)
          // coeff = dot(grad_Na, sigma @ grad_Nb)  # scalar
          // Kab = coeff * Ve * I3  # (3x3)

          // # Insert Kab into Kgeo_local at position [3a:3a+3, 3b:3b+3]
          // Kgeo_local[3*a:3*a+3, 3*b:3*b+3] += Kab

  
  
//}

//// K GEO - B MATRIX APPROACH
dev_t void Domain_d::CalcGeomStiffElementMatrix(){
  
  par_loop(e, m_elem_count){
   //# Insert Kab into 12x12 Kgeo_local at block [a][b]
  Matrix Bgeo(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?  

  Matrix Sigma(m_dim*2, m_dim *2);
  
// # Inputs:
// # - sigma: Cauchy stress tensor (3x3) [[σxx, σxy, σxz], [σxy, σyy, σyz], [σxz, σyz, σzz]]
// # - grad_N: List of shape function gradients [∇N1, ∇N2, ∇N3, ∇N4], each ∇Na is (3x1)
// # - Ve: Element volume

// # Step 1: Construct stress matrix Σ (6x6)
// Sigma = np.zeros((6, 6))
// # Fill diagonal blocks (3x3) of Σ
// Sigma[0:3, 0:3] = sigma[0, 0] * np.eye(3)  # σxx block
// Sigma[0:3, 3:6] = sigma[0, 1] * np.eye(3)  # σxy block
// Sigma[0:3, 6:9] = sigma[0, 2] * np.eye(3)  # σxz block
// # ... (repeat for σyy, σyz, σzz, ensuring symmetry)

  int offset_s = e;   //SCALAR
  int offset_t = offset_s * 6 ; //SYM TENSOR
  tensor3 sigma;
  sigma = FromFlatSym(m_sigma, offset_t );

// # Step 2: Build B_geo (6x12 matrix)
// B_geo = np.zeros((6, 12))
    
    for (int a=0;a<m_nodxelem;a++){
      double dNa[3];
      
    }
    // for a in range(4):
    // dN = grad_N[a]  # ∇Na (3x1)
    // B_geo[:, 3*a:3*a+3] = [
        // [dN[0], 0, 0],
        // [0, dN[1], 0],
        // [0, 0, dN[2]],
        // [dN[1], dN[0], 0],
        // [dN[2], 0, dN[0]],
        // [0, dN[2], dN[1]]
    // ]

// # Step 3: Compute K_geo (12x12)
// K_geo = Ve * B_geo.T @ Sigma @ B_geo

  }//ELEMENT 
}

////// K GEOMETRIC
///Kij_GEO = bTi sigma bj I3x3 dV
/// FOR EACH NODE PAIR
///dNdNiT sig dNdJ I3x3

//~ dev_t void Domain_d::assemblyForces(){

    //~ //if ()
  //~ //par_loop(n, m_node_count){
  //~ for (int n=0;n<m_node_count;n++){
    //~ for (int d=0;d<m_dim;d++)
      //~ m_fi[n*m_dim + d] = 0.0;
      
      //~ //printf("--------\n");    
      //~ for (int e=0; e<m_nodel_count[n];e++) {
        //~ int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
        //~ int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        //~ int offset  = eglob * m_nodxelem * m_dim;

        //~ for (int d=0;d<m_dim;d++){
          //~ //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
          //~ //if (n==9)
          //~ //  printf("%6e %6e %6e\n",m_fi[n*m_dim],m_f_elem[n*m_dim+1],m_f_elem[n*m_dim+2]);
          //~ m_fi[n*m_dim + d] += m_f_elem[offset + ne*m_dim + d];
        //~ }
          //~ if(m_thermal){
            //~ T[n] += dt * m_dTedt[eglob*m_nodxelem+ne];
	  //~ }
      //~ }
      //~ if (m_gp_count == 1 ) {  
        //~ for (int e=0; e<m_nodel_count[n];e++) {
          //~ int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
          //~ int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
          //~ int offset  = eglob * m_nodxelem * m_dim;
          //~ ////printf("glob %d, loc %d \n",n,ne);
          //~ for (int d=0;d<m_dim;d++){
            //~ //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
            //~ m_fi[n*m_dim + d] -= m_f_elem_hg [offset + ne*m_dim + d];
          //~ }
        //~ }      
      //~ }
      //~ // printf ("force %f %f %f\n",m_fi[m_dim*n],m_fi[m_dim*n+1],m_fi[m_dim*n+2]);
    //~ } // element


//~ }//assemblyForcesNonLock


};
