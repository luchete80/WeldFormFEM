/*************************************************************************/
/*  Thermal.C                                                     */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/

#include "Domain_d.h"
#include "Matrix_temp.h"
#include <stdio.h>

#include "tensor.cuh"

namespace MetFEM {
  
//dTdt = -invM * K * T
void Domain_d::ThermalCalcs(){
  par_loop(e, m_elem_count) {
    Matrix* Kt = new Matrix(m_nodxelem, m_nodxelem);

    for (int i = 0; i < m_nodxelem; ++i) {
      for (int j = 0; j < m_nodxelem; ++j) {
        double k = 0.0;
        for (int d = 0; d < m_dim; ++d)
          k += getDerivative(e, 0, d, i) * getDerivative(e, 0, d, j);
        Kt->Set(i, j, k*mat[e]->k_T*vol[e]);
      }
    }

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
    double heat = 0.9 * m_q_plheat[e];  // [J/m³]
    double total_heat = heat * vol[e];               // [J] total energy for element
    double heat_per_node = total_heat / m_nodxelem;  // even split [J/node]
        
    MatMul(*Kt, *Te, dTde);  // dTde = Kt * Te
    
    for (int i = 0; i < m_nodxelem; ++i)  m_dTedt[e*m_nodxelem + i] = 0.0;
    
    for (int i = 0; i < m_nodxelem; ++i) {
      int node_id = m_elnod[offset + i];
        double m_inv = 1.0 / m_mdiag[node_id];  // You must compute lumped thermal mass!
        // dTde = element heat flow vector
        // m_inv = 1/(lumped thermal mass) => m = rho * cp * V
        m_dTedt[e*m_nodxelem + i]  = - m_inv * dTde->getVal(i, 0);  // Include time step and minus sign

    //~ if (e==0){
    //~ //if (m_dTedt[e*m_nodxelem + i]>1.0e5 || std::isnan(m_dTedt[e*m_nodxelem + i]) ) {
        //~ printf("DEBUG - Element %d:\n", e);
        //~ printf("  m_q_plheat[e] = %.3e [J/m³]\n", m_q_plheat[e]);
        //~ printf("  vol[e] = %.3e [m³]\n", vol[e]);
        //~ printf("  heat_per_node = %.3e [J]\n", heat_per_node);
        //~ printf("  m_mdiag[node_id] = %.3e\n", m_mdiag[node_id]);
        //~ printf("  mat[e]->cp_T = %.3e\n", mat[e]->cp_T);
        //~ printf("  mat[e]->k_T = %.3e\n", mat[e]->k_T);
        //~ printf("  dT increment = %.3e [°C]\n", heat_per_node / (m_mdiag[node_id]*mat[e]->cp_T));
        //~ printf("  dTdTe = %.4e\n",m_dTedt[e*m_nodxelem + i]);
        //~ printf("  dTde: %.4e, minv: %.4e\n",dTde->getVal(i, 0),m_inv);
        
        //~ //exit(1);
    
    //~ }        
      
        m_dTedt[e*m_nodxelem + i] += heat_per_node / (m_mdiag[node_id]*mat[e]->cp_T);// (J) / (J/°C) = °C
        //printf ("temp inc  x node %.3e\n",heat_per_node / (m_mdiag[node_id]*mat[e]->cp_T));
    //~ if (m_dTedt[e*m_nodxelem + i]>1.0e10) { // Solo primer elemento



    }


		//~ f = 1./(d * Particles[i]->cp_T ); //[ºC m^3/J]
    //~ Particles[i]->dTdt = f * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source + Particles[i]->q_plheat * pl_work_heat_frac + Particles[i]->q_cont_conv);	
    // dT/dt = hA(T-Ts) *1/(mc) 

    delete Kt;
    delete Te;
    delete dTde;
  }//FOR ELEMENT
  
  
  ////// ATTENTION! ASSUMING CONSTANT CP, TO MODIFY
  par_loop(n,m_node_count){

    T[n] += q_cont_conv[n]*1.0/(m_mdiag[n] * mat[0]->cp_T)*dt;
    //printf("Node %d T %f , cond %f\n", n, T[n],q_cont_conv[n]);    
  }

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
