#ifndef _DOMAIN_H_
#define _DOMAIN_H_

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
  
  //// ELEMENT DATA
  double *p;
  double *str_rate;
  double *str_inc;
  
  
 __device__ inline void UpdateDensity(double dt);
 
 void AddBoxLength(const int &tag, double3 V, double Lx, double Ly, double Lz, double &r);
 void AllocateNodes(const int node_count);
};



#endif