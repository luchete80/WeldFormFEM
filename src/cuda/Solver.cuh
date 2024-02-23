#ifndef _SOLVER_CUH_
#define _SOLVER_CUH_

#include <cuda.h>
#include "cudautils.cuh"

#include "Domain_d.cuh"

namespace MetFEM{

class Solver {
public:
	Solver(){}
	
	virtual __host__ void Solve();
	virtual __host__ void Solve(Domain_d *dom){}
	
protected:
	Domain_d 	*m_dom;
};

class SolverChungHulbert:
public Solver {
public:
SolverChungHulbert(){}
SolverChungHulbert(Domain_d *d);

	virtual __host__ void Solve();

protected:


};

}; //Namespace

#endif