/************************************************************************

	Copyright 2012-2013 Luciano Buglioni

	Contact: luciano.buglioni@gmail.com

	This file is a part of FluxSol

	FluxSol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    FluxSol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#ifndef _SOLVER_H
#define _SOLVER_H

//#include "./Type/Vec3d.h"


//#include "Utils.h"

#include <time.h>


namespace MetFEM{

class Domain_d;

class Solver{
public:

	//Solver<number>():matdim(0)
  Solver():m_dof(0)
	{
	    maxiter=100;
      ptol=1.0e-3;
      //vtol=1.e-03;    //Todos los valores iguales
	}
	//Solver<number>(const int &d):
  Solver(const int &d):
	m_dof(d)
	{}
  void setDomain(Domain_d *d){m_dom = d;}  
  virtual void Allocate(){}
  virtual int Solve(){}
  virtual void SetRDOF(const int &, const double &val){}
  
  virtual void assemblyGlobalMatrix(){}
  
  virtual void Allocate(const int &dim){
    m_dof = dim;
  }
  virtual void applyDirichletBCs(){}

  virtual ~Solver(){}
  virtual inline const double & getU(int node, int dim) const {}

  virtual inline void addToU(int node, int dim, double delta) {};
protected:
  //Vec3D vtol;     //
  double ptol;    //
  double maxiter; //
  double actualiter;  //


	//number rtol,abstol;
  double rtol,abstol;
	int m_dof;		//CONST?
	int maxits;
  
  Domain_d *m_dom;

};

	//~ template <typename T>
    //~ void Solve(EqnSystem <T> &);

	//~ template <typename T>
    //~ void Solve(EqnSystem <T> &);
    
}

#endif
