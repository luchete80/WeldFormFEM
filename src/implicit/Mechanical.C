//CALCULATE INTERNAL FORCES
// F INT BT sigma

#include "Domain_d.h"
#include <iostream>
#include <vector>

// #include "tensor.cu"
#include "Matrix.h"

#if CUDA_BUILD
#include "tensor.cuh"
#else
//#include "Matrix.h"
#endif

#include "Tensor3.C"

using namespace std;

namespace MetFEM {

// 
// !!!!!! IT ASSUMES PRESSURE AND STRAIN RATES ARE ALREADY CALCULATED
// !!!!!! (AT t+1/2 to avoid stress at rigid rotations, see Benson 1992)
dev_t void Domain_d::CalcStressStrain(double dt){


  par_loop(e,m_elem_count){
      //printf("calculating sigma \n");
        // Jaumann rate terms
      tensor3 RotationRateT,SRT,RS;
      tensor3 RotRate;
      tensor3 StrRate;
      tensor3 ShearStress;
      tensor3 Sigma;

    //printf("calculating sigma %d\n", e);
    for (int gp=0;gp<m_gp_count;gp++){
      int offset_s = e * m_gp_count + gp;   //SCALAR
      int offset_t = offset_s * 6 ; //SYM TENSOR
      ShearStress = FromFlatSym(m_tau,          offset_t );
      StrRate     = FromFlatSym(m_str_rate,     offset_t );
      RotRate     = FromFlatAntiSym(m_rot_rate, offset_t );


      SRT = ShearStress * Trans(RotRate);
      RS  = RotRate * ShearStress;

      tensor3 test = StrRate-1.0/3.0*(StrRate.xx+StrRate.yy+StrRate.zz)*Identity();
      ShearStress	= ShearStress  + dt*(2.0* mat[e]->Elastic().G()*(StrRate - 1.0/3.0*Trace(StrRate) * Identity() ) + SRT+RS);
      
      double J2 = 0.5*(ShearStress.xx*ShearStress.xx +  2.0*ShearStress.xy*ShearStress.xy + 
                                      2.0*ShearStress.xz*ShearStress.xz + 
                     ShearStress.yy*ShearStress.yy+  
                                      2.0*ShearStress.yz*ShearStress.yz +               
                     ShearStress.zz*ShearStress.zz                 
                                     
                    );
      double sig_trial = sqrt(3.0*J2);

      if (sigma_y[e]<sig_trial){

        ShearStress = ShearStress * (sigma_y[e] / sig_trial);
        pl_strain[e] += (sig_trial - sigma_y[e]) / (3.0 *  mat[e]->Elastic().G());


      }

      Sigma = -p[offset_s] * Identity() + ShearStress;
 
      double Ep = 0;
			double dep=( sig_trial - sigma_y[e])/ (3.*mat[e]->Elastic().G() + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0

      ///// OUTPUT TO Flatten arrays
      ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
      //ToFlatSymPtr(Strain, 	strain,6*i);		
      ToFlatSymPtr(ShearStress, m_tau, offset_t);
      
    }//gp
  }//el < elcount


}

//To calculate before elemet Jacobian calc
dev_t void Domain_d::CalcElemVol(){
  par_loop(e,m_elem_count){
    double w;
    //TODO: CHANGE WEIGHT TO ARRAY
    if (m_gp_count == 1) {
      if (m_dim == 2) w = 4;//w = pow(2.0, m_dim);
      if (m_dim == 3)     
        if      (m_nodxelem == 4)  w = 1.0/6.0;
        else if (m_nodxelem == 8)  w = 8.0;
    } else                  w = 1.0;
    
    int offset = m_gp_count * e;
    vol[e] = 0.0;
    for (int gp=0;gp<m_gp_count;gp++){
      vol[e] += m_detJ[offset] * w;
    }  
    //if (e<10)
    //printf("Element %d Vol %f, det %f\n",e,vol[e],m_detJ[offset]);  
  }//el
}

////////////////////////////////
///// STIFNESS HOURGLASSING /////////////
////////////////////////////////
//~ Mode	Node 0	Node 1	Node 2	Node 3
//~ 1	1	-1	 0	 0
//~ 2	0	 1	-1	 0
//~ 3	0	 0	 1	-1
dev_t void Domain_d::calcElemHourglassForces()
{
  if (m_dim != 3 || m_nodxelem != 4) // Only 3D tetra mesh enters here
    return;

  int jmax = 3;  // 3 hourglass modes for tetra

  // Define hourglass mode vectors Sig(jmax x m_nodxelem)
  double sig_[3][4] = {
    {  1.0, -1.0,  0.0,  0.0 },
    {  0.0,  1.0, -1.0,  0.0 },
    {  0.0,  0.0,  1.0, -1.0 }
  };

  Matrix Sig(jmax, m_nodxelem);
  for (int i = 0; i < jmax; i++)
    for (int n = 0; n < m_nodxelem; n++)
      Sig.Set(i, n, sig_[i][n]);

  double hmod[3][3] = {0};  // m_dim x jmax

  par_loop(e, m_elem_count) {
    if (m_gp_count == 1) {
      int offset = e * m_nodxelem * m_dim;

      // Zero hourglass mode projections and forces
      for (int d = 0; d < m_dim; d++)
        for (int j = 0; j < jmax; j++)
          hmod[d][j] = 0.0;

      for (int d = 0; d < m_dim; d++)
        for (int n = 0; n < m_nodxelem; n++)
          m_f_elem_hg[offset + n * m_dim + d] = 0.0;


      for (int d = 0; d < m_dim; d++)
        for (int j = 0; j < jmax; j++)
          for (int n = 0; n < m_nodxelem; n++)
            hmod[d][j] += getVElem(e, n, d) * Sig.getVal(j, n);

      // Compute hourglass forces: stiffness-type resistance
      for (int d = 0; d < m_dim; d++)
        for (int j = 0; j < jmax; j++)
          for (int n = 0; n < m_nodxelem; n++)
            m_f_elem_hg[offset + n * m_dim + d] -= hmod[d][j] * Sig.getVal(j, n);

      // Hourglass stiffness coefficient (empirical tuning)
      double c_h = 0.1 * mat[e]->Elastic().E() * pow(vol[e], 2.0 / 3.0);  // E = Young's modulus

      // Scale HG forces
      for (int n = 0; n < m_nodxelem; n++)
        for (int d = 0; d < m_dim; d++)
          m_f_elem_hg[offset + n * m_dim + d] *= c_h;
    }
  }
}

void Domain_d::CalcElemIntForces(){
  par_loop(e, m_elem_count){
    //if (m_dim==3){
    Matrix B(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?
    B = getElemBMatrix(e);
    Matrix S(6,1);
    MatMul(MatMul(B.getTranspose(),S),B, m_Kmat[e]);
    
  }
}
};
