#include "Domain_d.h"

Domain_d::ThermalCalcs(){

//parallel loop here
for (int e=0)
  Matrix *Kt   = new Matrix(m_nodxelem,m_nodxelem );
  /////Tn+1=Tn+dt*invM(Q-KTn)
//BT x B
//faster than matmult? 

for (int i=0;i<m_nodxelem;i++)
  for (int j=0;j<m_nodxelem;j++){
    double k =0.0;
    for (int d=0;d<m_dim;d++)
      k +=getDerivative(e,0,i,d)*getDerivative(e,0,j,d);
    Kt->setVal(i,j,k);
  }
 }

Matrix *Te   = new Matrix(3,1);
Matrix *dTde   = new Matrix(3,1);
int offset=m_nodxelem*e;   ///// TO DO: improving var offset
for (int i=0;i<m_nodxelem;i++)
  Te[i] = T[m_elnod[offset+i]];
MatMul(*Kt, *Te, *dTde);


for (int i=0;i<m_nodxelem;i++){
   m_dTedt[offset+i]+=dTde->getVal(i); /////TODO: make function for this
 
}
  delete Kt, Te,dTde; 
} // for elements

}

Domain_d::calcInelasticHeatFraction(){

//parallel loop here
for (int e=0)
int offset_t = e * 6 ; //SYM TENSOR
tensor3 sig= FromFlatSym(m_sigma,          offset_t );


if (delta_pl_strain > 0.0) {
		tensor3 depdt = 1./dt*Strain_pl_incr;
		// Double inner product, Fraser 3-106
		//cout <<"depdt"<<endl;
		//cout << depdt<<endl;
		q_plheat 	= 
						m_sigma(0,0)*depdt(0,0) + 
						2.0*Sigma(0,1)*depdt(1,0) + 2.0*Sigma(0,2)*depdt(2,0) + 
						Sigma(1,1)*depdt(1,1) +
						2.0*Sigma(1,2)*depdt(2,1) + 
						Sigma(2,2)*depdt(2,2)
						;
		//cout << "plastic heat "<<q_plheat<<endl;
	}
  ps_energy[e] += q_plheat * dt;
}
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
