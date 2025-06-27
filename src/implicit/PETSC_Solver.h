/************************************************************************

	Copyright 2012-2015 Luciano Buglioni, Pablo Zitelli

	Contact:

	luciano.buglioni@gmail.com
	pzitelli@gmail.com

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

#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

#include "Solver.h"
#include "petscksp.h"
#include <vector>

namespace FluxSol

{

//TO MODIFY: MUST INCLUDE A EQN SYSTEM
//THIS ONTE CAN BE TEMPLATIZeD WITH MAtrIX; SPARSE MATRIX, MAT
//INHERIT FROM THIS SOLVER, SEQUENTIAL AND PARALLEL SOLVERS
template <typename number>
class PETSC_KSP_Solver:
public Solver<number>
{

	protected:

		// PETSC variables
		Mat A;
		Vec b; // right hand side
		Vec x; // solution, residual vectors
		Mat SysMat;

        PetscBool permute;
		MatOrderingType mat_ord_type;       //For Matrix reordering (See PETSC example 1 and 18)
		IS rowperm,colperm;

		KSP ksp; // linear solver context
		PC pc; // preconditioner context

		PetscErrorCode ierr;

		PetscMPIInt size;
		PetscScalar neg_one = -1.0,one = 1.0;
		PetscBool nonzeroguess = PETSC_FALSE;

        MPI_Comm       comm;

        //TO MODIFY
		//CAN HAVE A DOFHANDLER
		//DoFHandler<dim> &dofhandler;

		//TO MODIFY
		//REFERENCE TO DOFHANDLER OR GRID!!

	public:

	void PETSC_Init();
	//Constructors
	PETSC_KSP_Solver<number>():Solver<number>(0){};
	PETSC_KSP_Solver<number>(const int &d);

	void PreAllocateRows(const vector<int> &nnz);
	void PreAllocateRows(const PetscInt &cols);

	Mat Matrix(){return this->A;}

	void Solve();
	template <typename T>
    void Solve(EqnSystem <T> &TEqn){};

	void InsertRow(const int &row, const std::vector<int> &cols, const std::vector <double> &vals);

	void SetPermute(const PetscBool & b) {permute=b;}


	inline void SetMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value);
	inline void AddMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value);
	void SetbValues(const PetscInt &row, const PetscScalar &value);	//TO MODIFY, TEMPLATIZE TO VALUE
	void SetbValues(const vector<int> &row, const vector<number> &val);
	void SetxValues(const vector<int> &row, const vector<number> &val);
	void AddbValues(const PetscInt &row, const PetscScalar &value);	//TO MODIFY, TEMPLATIZE TO VALUE
	void SetbValues(const PetscScalar &value);

	void SetxValues(const PetscInt &row, const PetscScalar &value);	//TO MODIFY, TEMPLATIZE TO VALUE


	const std::vector <number> X() const;	//Returns X Solution
	const number X(const int &i) const;	//Returns X Solution
	const std::vector <number> B() const;
	const number B(const int &i) const;

	void ViewInfo();

	void Flush(){
		ierr = MatAssemblyBegin(this->A,MAT_FLUSH_ASSEMBLY);
		ierr = MatAssemblyEnd  (this->A,MAT_FLUSH_ASSEMBLY);
	}

	inline void ClearMat();

	void Destroy()
	{
        //cout << "Destying  ksp, "<< KSPDestroy(&ksp);

        MatDestroy(&A);
        //MatDestroy(&SysMat);
        KSPDestroy(&ksp);
        //PCDestroy(&pc);
        VecDestroy(&b);
        VecDestroy(&x);
	}
	~PETSC_KSP_Solver<number>()
	{
	    //this->Destroy();
	}

	void ShowInfo()
	{
     PetscReal   norm,norm2;
     PetscViewer viewer;
     Vec         res;
     this->comm = PETSC_COMM_WORLD;
     this->ierr = PetscViewerASCIIOpen(comm, "rhs.m", &viewer);CHKERRQ(this->ierr);
     this->ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(this->ierr);
     this->ierr = VecView(this->b,viewer);CHKERRQ(this->ierr);
     this->ierr = PetscViewerDestroy(&viewer);
     this->ierr = VecNorm(this->b, NORM_2, &norm2);CHKERRQ(this->ierr);

     this->ierr = PetscViewerASCIIOpen(comm, "solution.m", &viewer);CHKERRQ(this->ierr);
     this->ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(this->ierr);
     this->ierr = VecView(this->x,viewer);CHKERRQ(this->ierr);
     this->ierr = PetscViewerDestroy(&viewer);

     this->ierr = VecDuplicate(this->x, &res);CHKERRQ(this->ierr);
     this->ierr = MatMult(this->A, this->x, res);CHKERRQ(this->ierr);
     this->ierr = VecAXPY(this->b, -1.0, res);CHKERRQ(this->ierr);
     this->ierr = VecDestroy(&res);CHKERRQ(this->ierr);
     this->ierr = VecNorm(this->b,NORM_2,&norm);CHKERRQ(this->ierr);
     PetscPrintf(PETSC_COMM_WORLD,"[%d]%s |b-Ax|/|b|=%e, |b|=%e\n",0,__FUNCT__,norm/norm2,norm2);

     this->ierr = PetscViewerASCIIOpen(comm, "residual.m", &viewer);CHKERRQ(this->ierr);
     this->ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(this->ierr);
     this->ierr = VecView(this->b,viewer);CHKERRQ(this->ierr);
     this->ierr = PetscViewerDestroy(&viewer);

	}

	// ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
	// ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
	// ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	// /*
	// Always call PetscFinalize() before exiting a program. This routine
	// - finalizes the PETSc libraries as well as MPI
	// - provides summary and diagnostic information if certain runtime
	// options are chosen (e.g., -log_summary).
	// */
	// ierr = PetscFinalize();

};//PETSC Solver

} //FluxSol
#endif
