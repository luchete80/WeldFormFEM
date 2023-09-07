#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "NodeData.cuh"
#include "ElementData.cuh"


//DEVICE 
class Domain{
  public:
  //// NODE DATA, SEPARATE?
  double3* x; //Vector is double
	double3* v;
	double3* a;
	double3* u;
  
  int dim, nodxelem; //TODO: CHANGE BY ELEMENT 
  int node_count,elem_count;
  
  unsigned long *elnod;
  
  ElementData elem;
  NodeData    nod;
  
  int redint; //!Reduced integration
  
 __device__ inline void UpdateDensity(double dt);
 
 void AddBoxLength(const int &tag, double3 V, double Lx, double Ly, double Lz, double &r);
 void AllocateNodes(const int &);
 void AllocateElements(const int &, const int &);
};



#endif