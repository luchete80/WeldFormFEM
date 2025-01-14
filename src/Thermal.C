#include "Domain_d.h"

Domain_d::ThermalCalcs(){

//parallel loop here
for (int e=0)
  Matrix *Kt   = new Matrix(3,3);
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
int offset=m_nodxelem*e;
for (int i=0;i<m_nodxelem;i++){
   m_dTedt[offset+i]=; /////TODO: make function for this
 
}
  delete Kt; 
} // for element
