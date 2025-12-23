/*************************************************************************/
/*  Thermal.C                                                    */
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
#include "Matrix_temp.h"
#include <stdio.h>

#include "tensor.cuh"

namespace MetFEM {
  
//dTdt = -invM * K * T
void Domain_d::ThermalCalcs(){
  double w;
  if (m_gp_count == 1) {
    if (m_dim == 2){
      if      (m_nodxelem == 4) w = 4;//w = pow(2.0, m_dim);
      else if (m_nodxelem == 3) w = 1.0/2.0;
    }
    if (m_dim == 3){ 
      if      (m_nodxelem == 4)  w = 1.0/6.0;
      else if (m_nodxelem == 8)  w = 8.0;
    }
  }

  par_loop(e, m_elem_count) {
    Matrix* Kt = new Matrix(m_nodxelem, m_nodxelem);

    for (int i = 0; i < m_nodxelem; ++i) {
      for (int j = 0; j < m_nodxelem; ++j) {
        double k = 0.0;
        for (int d = 0; d < m_dim; ++d){
          double dNdxi_detJ = getDerivative(e, 0, d, i); //getDerivative is dNdx . detJ
          double dNdxj_detJ = getDerivative(e, 0, d, j);//getDerivative is dNdx . detJ
          k += dNdxi_detJ * dNdxj_detJ; 
        }
        Kt->Set(i, j, k*mat[e]->k_T/m_detJ[e]*w);
      }//i
    }//j
    //[K] = W/(m.K) . (1/m2) . m3 = W/K
    //[rho c V ] = kg/m3 . J/(kg.K) . m3  = J/K

    Matrix* Te = new Matrix(m_nodxelem, 1);
    Matrix* dTde = new Matrix(m_nodxelem, 1);
    
    //printf("Calc T elem\n");
    int offset = m_nodxelem * e;
    for (int i = 0; i < m_nodxelem; ++i) {
      Te->Set(i, 0, T[m_elnod[offset + i]]);
    }
    //Kt->Print();

		///// f = 1./(d * Particles[i]->cp_T ); //[ºC m^3/J]
    ///// Particles[i]->dTdt = f * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source + Particles[i]->q_plheat * pl_work_heat_frac + Particles[i]->q_cont_conv);	
    //[m_q_plheat] = [sigma:epsdot]= [N/m2].[1/s] = J/(m3.s) = W/m3
    double heat = 0.9 * m_q_plheat[e];  // [W/m³s]
    double elem_pow = heat * vol[e];               // [W] total energy for element
    double pow_per_node = elem_pow / m_nodxelem;  // even split [J/node]
        
    MatMul(*Kt, *Te, dTde);  // dTde = Kt * Te
    
    for (int i = 0; i < m_nodxelem; ++i)  m_dTedt[e*m_nodxelem + i] = 0.0;
    
    //~ for (int i = 0; i < m_nodxelem; ++i) {
      //~ int node_id = m_elnod[offset + i];
        //~ double m_inv = 1.0 / (m_mdiag[node_id]*mat[e]->cp_T);  // You must compute lumped thermal mass!
        //~ // dTde = element heat flow vector
        //~ // m_inv = 1/(lumped thermal mass) => m = rho * cp * V
        //~ m_dTedt[e*m_nodxelem + i]  = m_inv * ( pow_per_node - dTde->getVal(i, 0) );  // Include time step and minus sign
    //~ }


    for (int i = 0; i < m_nodxelem; ++i) {
      int node_id = m_elnod[offset + i];
        double m_inv = 1.0 / m_mdiag[node_id];  // You must compute lumped thermal mass!
        m_dTedt[e*m_nodxelem + i]  = - m_inv * dTde->getVal(i, 0);  // Include time step and minus sign

        m_dTedt[e*m_nodxelem + i] += pow_per_node / (m_mdiag[node_id]*mat[e]->cp_T);// (J) / (J/°C) = °C

    }
    
    delete Kt;
    delete Te;
    delete dTde;
  }//FOR ELEMENT
  
  double dTdt[m_node_count];
  
  //CHECK dT!!
  // dt < rho c h^2 / (2k)
  //
  ////////// ASSEMBLY
  ////// ATTENTION! ASSUMING CONSTANT CP, TO MODIFY
  //double max_dTdt = 100.0;
  par_loop(n,m_node_count){
      dTdt[n] = 0;

      for (int e=0; e<m_nodel_count[n];e++){ 
        int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
        int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        int offset  = eglob * m_nodxelem * m_dim;
        
        dTdt[n] += m_dTedt[eglob*m_nodxelem+ne];
        
    }
    //if (dTdt[n]> max_dTdt) dTdt[n]= max_dTdt; 
    //if (dTdt[n]<-max_dTdt) dTdt[n]=-max_dTdt;
        
    T[n] += (dTdt[n] + q_cont_conv[n]*1.0/(m_mdiag[n] * mat[0]->cp_T) ) *dt;

  }//Node Loop

}

void Domain_d::calcInelasticHeatFraction(){

  par_loop(e, m_elem_count) {
  int offset_t = e * 6 ; //SYM TENSOR
  tensor3 sig             = FromFlatSym(m_sigma,          offset_t );
  tensor3 Strain_pl_incr  = FromFlatSym(m_strain_pl_incr, offset_t );

  //TODO: ADD A PLASTIC STRAIN INCREMENT FLAG
  //~ if (delta_pl_strain > 0.0) {
      
  tensor3 depdt = 1./dt * Strain_pl_incr;
      //~ // // Double inner product, Fraser 3-106
      //printf("depdt %.3e\n", depdt);
      //~ // //cout << depdt<<endl;
      m_q_plheat[e] 	= 
        sig.xx * depdt.xx + 
        2.0*sig.xy * depdt.yx + 2.0 * sig.xz*depdt.zx +              //~ // ;
        sig.yy*depdt.yy +      //~ // //cout << "plastic heat "<<q_plheat<<endl;
        2.0*sig.yz*depdt.yz +     //~ // ps_energy[e] += q_plheat * dt;
        sig.zz*depdt.zz; //Parallel
      //~ if (m_q_plheat[e]>1.0e-5)
        //~ printf("m_q_plheat %.3e\n",m_q_plheat[e]);
    }
    
  //~ #ifndef  CUDA_BUILD
  //~ for (int e=0;e<m_elem_count;e++)
    //~ m_pl_energy +=m_q_plheat[e]*dt;
  //~ #else

  //~ #endif
}//IHF

////// Dtotal = Dmech + Dthermal
///// This substract thermal to work Stress with Mechanical part
dev_t void Domain_d::calcThermalExpansion(){
  par_loop(e, m_elem_count) {
    int offset_t = 6 *e;

    double dTdt_gp = 0.0;
    for (int i = 0; i < m_nodxelem; ++i) {
        int node = m_elnod[e*m_nodxelem + i];
        dTdt_gp += 1.0/m_nodxelem * m_dTedt[node];
    }

    tensor3 StrRate     = FromFlatSym(m_str_rate,     offset_t );
    StrRate = StrRate - mat[e]->exp_T * dTdt_gp * Identity();
    ToFlatSymPtr(StrRate,m_str_rate,   offset_t);
  }
}


/*

inline void Particle::CalcPlasticWorkHeat(const double &dt){
	
	if (delta_pl_strain > 0.0) {
		Mat3_t depdt = 1./dt*Strain_pl_incr;
		// Double inner product, Fraser 3-106
		//cout <<"depdt"<<endl;
		//cout << depdt<<endl;
		q_plheat 	= 
						Sigma(0,0)*depdt(0,0) + 
						2.0*Sigma(0,1)*depdt(1,0) + 2.0*Sigma(0,2)*depdt(2,0) + 
						Sigma(1,1)*depdt(1,1) +
						2.0*Sigma(1,2)*depdt(2,1) + 
						Sigma(2,2)*depdt(2,2)
						;
		//cout << "plastic heat "<<q_plheat<<endl;
	}
  ps_energy += q_plheat * dt;
}

*/

};
