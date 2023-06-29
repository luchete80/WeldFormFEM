#ifndef _DOMAIN_H_
#define _DOMAIN_H_

class Domain{
  
  //// NODE DATA, SEPARATE?
  double3* x; //Vector is double
	double3* v;
	double3* a;
	double3* u;
 
 
 __device__ inline void UpdateDensity(double dt);
  
};



#endif