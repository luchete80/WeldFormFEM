#include "Domain_d.h"

Domain_d::ThermalCalcs(){
Matrix *Kt   = new Matrix(3,3);
/////Tn+1=Tn+dt*invM(Q-KTn)
//BT x B
//faster than matmult? 

for (int i=0;i<m_nodxelem;i++)
  for (int j=0;j<m_nodxelem;j++){
    double e =0.0;
    for (int d=0;d<m_dim;d++)
      e+=
  }
}
