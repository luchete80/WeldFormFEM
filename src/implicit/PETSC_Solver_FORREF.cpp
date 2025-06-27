/************************************************************************

Copyright 2012-2013 Luciano Buglioni

Contact:
    luciano.buglioni@gmail.com

This file is a part of FluxSol

For a copy of the GNU General Public License,
see <http://www.gnu.org/licenses/>.

*************************************************************************/

#include "PETSC_Solver.h"

namespace MetFEM
{

template <typename number>
void PETSC_KSP_Solver<number>::PETSC_Init()
{
	int n=this->matdim;

	cout << "Initializing Solver..."<<endl;

	//If no values passed
	int argc;
	char **args;

	char help[100];

    //NOT TO USE THIS!! IT IS CRASHING
	argc=0;
	args=NULL;


	cout << "Initializing PETSC"<<endl;
	PetscInitialize(&argc,&args,(char *)0,help);

    PetscMPIInt size,rank;

    cout << "PETSC Initialized"<<endl;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
	//MPI_Init()

	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-nonzero_guess",&nonzeroguess,PETSC_NULL);CHKERRQ(ierr);

	/*
	Create vectors. Note that we form 1 vector from scratch and
	then duplicate as needed.
	*/
	cout <<"Creating matrix and vectors"<<endl;
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	//ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

	/*
	Create matrix. When using MatCreate(), the matrix format can
	be specified at runtime.
	Performance tuning note: For problems of substantial size,
	preallocation of matrix memory is crucial for attaining good
	performance. See the matrix chapter of the users manual for details.
	*/
    //THIS FUNCTIONS MUST NOT BE CALLED (SLOW ACCORDING TO PETSC MANUAL)
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);	//Instead of create Sij
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);




	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	//TEMPORARY
	//THIS IS NOT RECOMMENDED DIRECTLY
	//ierr = MatCreateSeqBAIJ(PETSC_COMM_WORLD,PetscInt bs,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A)
	//SEQSBAIJ

	//Symmetric, THIS IS TEMP, TO MODIFY
	//ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);


    vector<int> diagonal_nonzeros, off_diagonal_nonzeros;
	int nextCellCount;

	int nVars=3;


//	// Calculate space necessary for matrix memory allocation
//	for (int cellit=grid[gid].cell.begin();cit!=grid[gid].cell.end();cit++) {
//		nextCellCount=0;
//		for (it=(*cit).faces.begin();it!=(*cit).faces.end();it++) {
//			if (grid[gid].face[*it].bc==INTERNAL_FACE) {
//				nextCellCount++;
//			}
//		}
//		for (int i=0;i<nVars;++i) {
//			diagonal_nonzeros.push_back( (nextCellCount+1)*nVars);
//			off_diagonal_nonzeros.push_back( ((*cit).ghosts.size())*nVars);
//		}
//	}
//
//	MatCreateMPIAIJ(
//					PETSC_COMM_WORLD,
//					grid[gid].cellCount*nVars,
//					grid[gid].cellCount*nVars,
//					grid[gid].globalCellCount*nVars,
//					grid[gid].globalCellCount*nVars,
//					0,&diagonal_nonzeros[0],
//					0,&off_diagonal_nonzeros[0],
//					&impOP);

    //  MatCreateSeqAIJ(MPI_Comm comm,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A)
//    	comm 	- MPI communicator, set to PETSC_COMM_SELF
//	m 	- number of rows
//	n 	- number of columns
//	nz 	- number of nonzeros per row (same for all rows)
//	nnz 	- array containing the number of nonzeros in the various rows (possibly different for each row) or NULL

	//Must call MatXXXSetPreallocation() or MatSetUp() on argument 1 "mat" before MatSetValues()!

	ierr = MatSetUp(A);
	//
	//ierr= MatSeqAIJSetPreallocation(Mat B,PetscInt nz,const PetscInt nnz[])
	//where
	//B	- The matrix-free
	//nz	- number of nonzeros per row (same for all rows)
	//nnz	- array containing the number of nonzeros in the various rows (possibly different for each row) or NULL


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create the linear solver and set various options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	Create linear solver context
	*/
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	/*
	Set operators. Here the matrix that defines the linear system
	also serves as the preconditioning matrix.
	*/

	ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	//ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	//KSPSetFromOptions(ksp);

	/*
	Set linear solver defaults for this problem (optional).
	- By extracting the KSP and PC contexts from the KSP context,
	we can then directly call any KSP and PC routines to set
	various options.
	- The following four statements are optional; all of these
	parameters could alternatively be specified at runtime via
	KSPSetFromOptions();
	*/
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	//ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
	//ierr = PCSetType(pc,PCGAMG);CHKERRQ(ierr);
	//ierr = PCSetType(pc,PCICC);CHKERRQ(ierr);
	//ierr = PCSetType(pc,PCICC);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
	//ierr = KSPSetTolerances(ksp,1.e-3,1.e-2,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	//PetscErrorCode  KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)

    //ksp	- the Krylov subspace context
    //rtol	- the relative convergence tolerance, relative decrease in the (possibly preconditioned) residual norm
    //abstol	- the absolute convergence tolerance absolute size of the (possibly preconditioned) residual norm
    //dtol	- the divergence tolerance, amount (possibly preconditioned) residual norm can increase before KSPConvergedDefault() concludes that the method is diverging
    //maxits	- maximum number of iterations to use

	KSPSetInitialGuessKnoll(ksp,PETSC_TRUE);
	//KSPSetType(ksp,KSPGMRES);
	KSPSetType(ksp,KSPBCGS);    //BiCGSTAB

	/*
	Set runtime options, e.g.,
	-ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	These options will override those specified above as long as
	KSPSetFromOptions() is called _after_ any other customization
	routines.
	*/
	//ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr); THIS IS PREFERRED BEFORE PREVIOUS SOLVER
	//AND MUST BE CALLED SETUP,
	if (nonzeroguess)
	{
		PetscScalar p = .5;
		ierr = VecSet(x,p);CHKERRQ(ierr);
		ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
	}


	ierr=VecSet(this->b,0.);
	ierr=VecSet(this->x,0.);

	//REORDErING VARS
	http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex18.c.html
    permute=PETSC_FALSE;     //Matrix reordering
    mat_ord_type=MATORDERINGRCM;
    rowperm = NULL,colperm = NULL;

	cout << "[I] Solver Initialized." <<endl;


}
//
template <typename number>
PETSC_KSP_Solver<number>::
PETSC_KSP_Solver(const int &d):
Solver<number>(d)
{
	PETSC_Init();
}


template <typename number>
void PETSC_KSP_Solver<number>::PreAllocateRows(const vector <int> &nnz)
{
	MatSeqAIJSetPreallocation(this->A,PETSC_NULL,&nnz[0]);
}

template <typename number>
void PETSC_KSP_Solver<number>::PreAllocateRows(const PetscInt &cols)
{
	MatSeqAIJSetPreallocation(this->A,cols,PETSC_NULL);
    //PetscErrorCode  MatSeqAIJSetPreallocation(Mat B,PetscInt nz,const PetscInt nnz[])
    //
    //B	- The matrix
    //nz	- number of nonzeros per row (same for all rows)
    //nnz	- array containing the number of nonzeros in the various rows (possibly different for each row) or NULL

}

template <typename number>
void PETSC_KSP_Solver<number>::Solve()
{

	ierr = MatAssemblyBegin(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(this->x);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(this->x);CHKERRQ(ierr);

    if (permute)
    {
        cout << "performing matrix permutation..."<<endl;
        Mat Aperm;
        MatGetOrdering(this->A,mat_ord_type,&rowperm,&colperm);
        MatPermute(this->A,rowperm,colperm,&Aperm);
        VecPermute(b,colperm,PETSC_FALSE);
        MatDestroy(&A);
        this->A    = Aperm;               /* Replace original operator with permuted version */
        //MatDestroy(&Aperm);
    }

	//ierr = MatSetOption(this->A,MAT_SYMMETRIC,PETSC_TRUE);


	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	/*
	View solver info; we could instead use the option -ksp_view to
	print this info to the screen at the coknclusion of KSPSolve().
	*/
	//ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Check solution and clean up
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	Check the error
	*/
	PetscReal norm;
	PetscInt its;

	ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %A, Iterations %D\n",
	norm,its);CHKERRQ(ierr);

	if (permute) {VecPermute(x,rowperm,PETSC_TRUE);}

	/*
	Free work space. All PETSc objects should be destroyed when they
	are no longer needed.
	*/

	ISDestroy(&rowperm);  ISDestroy(&colperm);
}

template <typename number>
void PETSC_KSP_Solver<number>::InsertRow(const int &row, const std::vector<int> &cols, const std::vector <double> &vals)
{
	//ierr = MatSetValues(this->A,1,row,cols.size(),cols,&vals[0],INSERT_VALUES);CHKERRQ(ierr);
}

template <typename number>
void PETSC_KSP_Solver<number>::SetMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value)
{
	ierr=MatSetValues(this->A,1,&row,1,&col,&value,INSERT_VALUES);
}

template <typename number>
void PETSC_KSP_Solver<number>::AddMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value)
{
	ierr=MatSetValues(this->A,1,&row,1,&col,&value,ADD_VALUES);
}


template <typename number>
void PETSC_KSP_Solver<number>::ViewInfo()
{

	ierr = MatAssemblyBegin(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	cout << "Matrix A: "<<endl;
	MatView(this->A,PETSC_VIEWER_STDOUT_SELF);
	cout << "RHS: "<<endl;
	VecView(b,PETSC_VIEWER_STDOUT_SELF);
	cout << "Solution: "<<endl;
	VecView(x,PETSC_VIEWER_STDOUT_SELF);
}

template <typename number>
void PETSC_KSP_Solver<number>::SetbValues(const PetscInt &row, const PetscScalar &value)
{
	ierr=VecSetValues(this->b,1,&row,&value,INSERT_VALUES);
}

template <typename number>
void PETSC_KSP_Solver<number>::SetxValues(const PetscInt &row, const PetscScalar &value)
{
	ierr=VecSetValues(this->x,1,&row,&value,INSERT_VALUES);
}

template <typename number>
void PETSC_KSP_Solver<number>::SetbValues(const vector<int> &row, const vector<number> &val)
{
	ierr=VecSetValues(this->b,1,&row[0],&val[0],INSERT_VALUES);
}

template <typename number>
void PETSC_KSP_Solver<number>::SetxValues(const vector<int> &row, const vector<number> &val)
{
	ierr=VecSetValues(this->x,1,&row[0],&val[0],INSERT_VALUES);
}


//template <typename number,int dim>
//void PETSC_KSP_Solver<number,dim>::AddbValues(const PetscInt &row, const PetscScalar &value)
//{
//	ierr=VecSetValues(this->b,1,&row,&value,ADD_VALUES);
//}
//
//template <typename number,int dim>
//void PETSC_KSP_Solver<number,dim>::SetbValues(const PetscScalar &value)
//{
//	ierr=VecSet(this->b,value);
//}
//



template <typename number>
const std::vector <number>
PETSC_KSP_Solver<number>::X() const
{
	std::vector <number> v(this->matdim);

	int ix[1];
	double y[1];

	for (int i=0;i<this->matdim;i++) //TO EVALUATE SPEED
	{
		ix[0]=i;
		VecGetValues(x,1,ix,y);
		number val;
		val=y[0];
		v[i]=val;
	}

	//PetscErrorCode  VecGetValues(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])

	return v;
}

template <typename number>
void PETSC_KSP_Solver<number>::ClearMat()
{
	//TO MODIFY
	int row,col;
	//for (int i=0;i<)
	//ierr=MatSetValue(this->A,1,&row,1,&col,&value,INSERT_VALUES);

}

template <typename number>
const number
PETSC_KSP_Solver<number>::B(const int &i)const
{
	double y[1];
	int ix[1];
	ix[0]=i;

	VecGetValues(this->b,1,ix,y);
	//cout << "value get" << y[0] <<endl;
	number n=y[0];
	return n;
}

template <typename number>
const number
PETSC_KSP_Solver<number>::X(const int &i)const
{
	double y[1];
	int ix[1];
	ix[0]=i;

	VecGetValues(this->x,1,ix,y);
	//cout << "value get" << y[0] <<endl;
	number n=y[0];
	return n;
}

template <typename number>
const std::vector <number>
PETSC_KSP_Solver<number>::B() const
{
	std::vector <number> v(this->matdim);

	int ix[1];
	double y[1];

	for (int i=0;i<this->matdim;i++) //TO EVALUATE SPEED
	{
		ix[0]=i;
		VecGetValues(this->b,1,ix,y);
		number val;
		val=y[0];
		v[i]=val;
	}

	//PetscErrorCode  VecGetValues(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])

	return v;
}


template <typename T>
void Solve(EqnSystem <T> &TEqn)
{


    clock_t ittime_begin, ittime_end,ittime_start;
    double ittime_spent;

    ittime_begin = clock();
    ittime_start = clock();

    ///// PETSC SOLVE MANUALLY ///////

        int numberofcomp=pow(3.,TEqn.Dim());
        int totrows=numberofcomp*TEqn.Num_Eqn();

    cout << "Creating Solver"<<endl;
        PETSC_KSP_Solver<double> Solver(totrows);

    cout << "Number of comps" << numberofcomp<<endl;
    cout << "Number of rows" << totrows<<endl;

    vector <int> nonzerosperrow;
        for (int e=0;e<TEqn.Num_Eqn();e++)      //Aca voy con las filas de a 2
    {
        //TO MODIFY
        for (int dim=0;dim<numberofcomp;dim++)  nonzerosperrow.push_back(5);
    }

    cout << "Allocating Rows..."<<endl;
    Solver.PreAllocateRows(nonzerosperrow);

    ittime_spent = (double)(clock() - ittime_begin) / CLOCKS_PER_SEC;

    cout << "PETSC Creating Solver and Preallocating time: "<<ittime_spent<<endl;

    ittime_begin = clock();


    cout << "Assembying Eqns"<<endl;
        for (int e=0;e<TEqn.Num_Eqn();e++)      //Aca voy con las filas de a 2
        {
            //Width Assign

        //cout << "Assemblying Eqn "<<e<<endl;
                vector <double> ap=TEqn.Eqn(e).Ap().Comp();
                Scalar ap_sc=TEqn.Eqn(e).Ap();
                Scalar value;

                int row=e*numberofcomp;

        vector <double> col;
        //CENTRAL COEFFS
        col=ap;
        for (int dim=0;dim<numberofcomp;dim++)
            Solver.SetMatVal(row+dim, row+dim, col[0]);    //An is scalar


                //Look trough entire width for neighbours id
                //The main idea is to look through eqn width
                //Real cell id is taken, and then are watched all neighbours to check if each real cell id belongs to neighbours vector
                for (int nc=0;nc<TEqn.Eqn(e).Num_Neighbours();nc++)
                {
                        int realcellid=TEqn.Eqn(e).Neighbour(nc);   //Wich cell

                        col=TEqn.Eqn(e).An(nc).Comp();
                        int columnid=numberofcomp*realcellid;

            //cout << "Found Cell " <<endl;
            for (int dim=0;dim<numberofcomp;dim++)
                Solver.SetMatVal(row+dim, columnid+dim, col[0]);    //An is scalar

                }//En of neighbours


        }//End of cells for (int e=0;e<TEqn.Num_Eqn();e++)      //Aca voy con las filas de a 2




    //cout << "R vector (from zero)"<<endl;
        //V_SetAllCmp(&R,0.0);
        for (int e=0;e<TEqn.Num_Eqn();e++)
        {
            //cout << "Eqn " << e<<endl;
            //cout << "[" <<e<<"]: "  ;
                vector <double> source=TEqn.Eqn(e).Source().Comp();
                for (int dim=0;dim<numberofcomp;dim++)
        {
            Solver.SetbValues(e*numberofcomp+dim, source[dim]);
        }
        //cout << endl;
        }


    //cout << "Assemblying vector"<<endl;
    //Initial values
    //vector <double> r(TEqn.Num_Eqn()*numberofcomp);
//    double val;
//    int row;
//         for (int e=0;e<TEqn.Num_Eqn();e++)
//         {
//         //cout << "e= "<<e<<endl;
//         for (int dim=0;dim<numberofcomp;dim++)
//         {
//             row=numberofcomp*e+dim;
//             //val=0.5*TEqn.InitField().Val(e).Comp()[dim];
//             val=0.;
//             Solver.SetxValues(row,val);
//             //r[numberofcomp*e+dim]=TEqn.InitField().Val().Comp[dim];
//             //r[dim]=Solver.SetxValues(numberofcomp*e+dim);
//             //cout << "xi= "<<numberofcomp*e+dim<<", ";
//             //cout <<U.Cmp[numberofcomp*e+dim+1]<<" ";
//         }
//       // TEqn.Eqn(e).X()=r;
//
//     }

     cout << "Assembled"<<endl;

    ittime_spent = (double)(clock() - ittime_begin) / CLOCKS_PER_SEC;
        cout << "Assemblying elapsed time: "<<ittime_spent<<endl;

    ittime_begin = clock();

    Solver.Solve();

    ittime_spent = (double)(clock() - ittime_begin) / CLOCKS_PER_SEC;

    cout << "PETSC Solving elapsed time: "<<ittime_spent<<endl;

    ittime_begin = clock();
    cout <<"Solver Results "<<endl;
    vector <double> r(numberofcomp);
         for (int e=0;e<TEqn.Num_Eqn();e++)
         {
         //cout << "e= "<<e<<endl;
         for (int dim=0;dim<numberofcomp;dim++)
         {
             r[dim]=Solver.X(numberofcomp*e+dim);
             //cout << "xi= "<<numberofcomp*e+dim<<", ";
             //cout <<U.Cmp[numberofcomp*e+dim+1]<<" ";
         }
        TEqn.Eqn(e).X()=r;

     }

    ittime_spent = (double)(clock() - ittime_begin) / CLOCKS_PER_SEC;
        cout << "PETSC Allocating results time: "<<ittime_spent<<endl;

    ittime_begin = clock();

    //Solver.ShowInfo();


     cout << "Destroying "<<endl;

    Solver.Destroy();

    ittime_spent = (double)(clock() - ittime_begin) / CLOCKS_PER_SEC;
        cout << "PETSC Destroy time: "<<ittime_spent<<endl;

    cout << "Destroyed"<<endl;


    ittime_spent = (double)(clock() - ittime_start) / CLOCKS_PER_SEC;
        cout << "PETSC Total time: "<<ittime_spent<<endl;

}

} //FluxSol

#include "PETSC_Solver.inst"
