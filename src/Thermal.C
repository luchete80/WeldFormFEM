#include "Domain_d.h"
#include "Matrix_temp.h"
#include <stdio.h>

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
    
    MatMul(*Kt, *Te, dTde);  // dTde = Kt * Te

    for (int i = 0; i < m_nodxelem; ++i) {
      int node_id = m_elnod[offset + i];
        double m_inv = 1.0 / m_mdiag[node_id];  // You must compute lumped thermal mass!
        // dTde = element heat flow vector
        // m_inv = 1/(lumped thermal mass) => m = rho * cp * V
        m_dTedt[e*m_nodxelem + i] += -dt * m_inv * dTde->getVal(i, 0);  // Include time step and minus sign
    }


		//~ f = 1./(d * Particles[i]->cp_T ); //[ºC m^3/J]
    //~ Particles[i]->dTdt = f * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source + Particles[i]->q_plheat * pl_work_heat_frac + Particles[i]->q_cont_conv);	
    // dT/dt = hA(T-Ts) *1/(mc) 

    delete Kt;
    delete Te;
    delete dTde;
  }//FOR ELEMENT
  
  
  par_loop(n,m_node_count){

    T[n] += q_cont_conv[n]*1.0/(m_mdiag[n] * 960.0)*dt;
    //printf("Node %d T %f , cond %f\n", n, T[n],q_cont_conv[n]);    
  }

}

void Domain_d::calcInelasticHeatFraction(){

  // //parallel loop here
  // for (int e=0)
  // int offset_t = e * 6 ; //SYM TENSOR
  // tensor3 sig= FromFlatSym(m_sigma,          offset_t );


  //~ if (delta_pl_strain > 0.0) {
      
      //~ tensor3 depdt = 1./dt*Strain_pl_incr;
      //~ // // Double inner product, Fraser 3-106
      //~ // //cout <<"depdt"<<endl;
      //~ // //cout << depdt<<endl;
      //~ // q_plheat 	= 
              //~ // m_sigma(0,0)*depdt(0,0) + 
              //~ // 2.0*Sigma(0,1)*depdt(1,0) + 2.0*Sigma(0,2)*depdt(2,0) + 
              //~ // Sigma(1,1)*depdt(1,1) +
              //~ // 2.0*Sigma(1,2)*depdt(2,1) + 
              //~ // Sigma(2,2)*depdt(2,2)
              //~ // ;
      //~ // //cout << "plastic heat "<<q_plheat<<endl;
    //~ // }
    //~ // ps_energy[e] += q_plheat * dt;
   //~ }
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
