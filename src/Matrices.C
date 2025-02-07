#include "Domain_d.h"

#include "Matrix_temp.h"

#ifdef CUDA_BUILD
#include <cuda_runtime.h>
#else
#endif

namespace MetFEM {
  
// __global__ void histogram(int* color, int* buckets)
// {
 // int i = threadIdx.x + blockDim.x * blockIdx.x;
 // int c = colors[i];
 // atomicAdd(&buckets[c], 1);
// }

// TODO: REPLACE ATOMICADD PER ANOTHER METHOD


///////HERE THE LOOP IS PER NODE INSTEAD THAN PER ELEMENT
  dev_t void Domain_d::assemblyForces(){

    //if ()
  par_loop(n, m_node_count){

      for (int e=0; e<m_nodel_count[n];e++) {
        for (int d=0;d<m_dim;d++)
          m_fi[n*m_dim + d] = 0.0;
      }
          
      for (int e=0; e<m_nodel_count[n];e++) {
        int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
        int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        int offset  = eglob * m_nodxelem * m_dim;
        ////printf("glob %d, loc %d \n",n,ne);
        for (int d=0;d<m_dim;d++){
          //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
          m_fi[n*m_dim + d] += m_f_elem[offset + ne*m_dim + d];
        }
          if(m_thermal){
            T[n] += dt * m_dTedt[eglob*m_nodxelem+ne];
	  }
      }
      if (m_gp_count == 1 ) {  
        for (int e=0; e<m_nodel_count[n];e++) {
          int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
          int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
          int offset  = eglob * m_nodxelem * m_dim;
          ////printf("glob %d, loc %d \n",n,ne);
          for (int d=0;d<m_dim;d++){
            //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
            m_fi[n*m_dim + d] -= m_f_elem_hg [offset + ne*m_dim + d];
          }
        }      
      
      }
      // printf ("force %f %f %f\n",m_fi[m_dim*n],m_fi[m_dim*n+1],m_fi[m_dim*n+2]);
    } // element

  }//assemblyForcesNonLock


__global__ void assemblyForcesKernel(Domain_d *dom_d){
	dom_d->assemblyForces();
}


/*
 dev_t void Domain_d::assemblyForces(){

   //if ()
   for (int n=0; n<m_node_count;n++) {
     for (int d=0;d<m_dim;d++)
       m_fi[n*m_dim + d] = 0.0;
   }
  
 par_loop(e, m_elem_count){


        
     for (int e=0; e<m_elem_count;e++) {
       //printf("glob %d, loc %d \n",n,ne);
       for (int ne=0;ne<m_nodxelem;ne++){
         int nglob = m_elnod[e * m_nodxelem+ne];
         int offset  = e * m_nodxelem * m_dim + ne*m_dim;
         for (int d=0;d<m_dim;d++){

           //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
           m_fi[3*nglob + d ] +=  m_f_elem[offset + d]
                                         - m_f_elem_hg[offset + d];
         }
       }
     }

     //printf ("force %f %f %f\n",m_fi[m_dim*n],m_fi[m_dim*n+1],m_fi[m_dim*n+2]);
   }  //element

 }
  
  //assemblyForcesNonLock

*/

//// ORIGINAL WAS ASSEMBLED BY ELEMENT
//subroutine assemble_mass_matrix ()
// m_glob (:,:) = 0.0d0
  // do e = 1, elem_count
    // !print *, "elem ", e
    // do n1 =1, nodxelem
      // do n2=1, nodxelem
            // !print *, "elem ", e, "node ", n, " i j matm ",i, j, elem%matm (e,dim*(n-1)+i,dim*(n2-1)+j)            
            // iglob  = elem%elnod(e,n1) 
            // jglob  = elem%elnod(e,n2) 
            // ! print *, "iloc, jloc ",dim*(n-1)+i, dim*(n2-1)+j, "iglob, jglob", iglob,jglob
            // m_glob(iglob,jglob) = m_glob(iglob,jglob) + elem%matm (e,n1,n2)
      // end do !dof1
    // end do ! dof2 
  // end do ! e
// end subroutine
// dev_t void Domain_d::assemblyMassMatrix(Domain_d *dom_d)
// {
  // par_loop(n, m_node_count){  
    // offset = n*m_node_count;
    // for (int i=0; i<m_node_count;i++) m_mglob[offset + i] = 0.0;
      // for (int e=0; e<m_nodel_count[n];e++) {
        // int eglob   = m_nodel     [m_nodel_offset[n]+e];
        // int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        
      // }//elcount 
  // }//par loop
// }
  
  //CPU TEST (SINCE IS CALULATED ONLY ONCE
  /// THIS IS WRITING AT SAME MEM LOCATION
  // ONLY FOR CPU, NOT CUDA
  // USES ELNOD_H; CHANGE TO GPU
  
  /// NOT WORKING
  void Domain_d::calcMassDiagFromElementNodes(const double &rho) // To use existing array
  {
    Matrix m_glob(m_node_count,m_node_count);
    Matrix temph(m_dim, m_nodxelem);
    Matrix tempm(m_nodxelem, m_nodxelem);
    int n1,n2, iglob,jglob;
    
    // m_glob (:,:) = 0.0d0
    for (int e=0;e<m_elem_count;e++){
      // !print *, "elem ", e


    int offset = e * m_nodxelem * m_nodxelem;
  
  
    //if (m_gp_count == 1 ) {  

      //double w = pow(2.0, m_dim);
      double f = 1.0/(pow(2.0,m_dim)); //THIS INCLUDES WEIGHT!!!
      for (int k=0;k<m_nodxelem;k++){  
        for (int d=0;d<m_dim;d++)
          temph.Set(d,k,f);
      }
      // TODO: DO NOT PERFORM THIS MULT
      
      MatMul(temph.getTranspose(),temph,&tempm);


      for( n1 =0; n1<m_nodxelem;n1++){ //LOCAL n1
        for( n2 =0; n2<m_nodxelem;n2++){ //LOCAL n2
              // !print *, "elem ", e, "node ", n, " i j matm ",i, j, elem%matm (e,dim*(n-1)+i,dim*(n2-1)+j)            
              int iglob  = elnod_h[e*m_nodxelem+n1];
              int jglob  = elnod_h[e*m_nodxelem+n2];              
              // jglob  = elnod(e,n2) 
              // ! print *, "iloc, jloc ",dim*(n-1)+i, dim*(n2-1)+j, "iglob, jglob", iglob,jglob
              // m_glob(iglob,jglob) = m_glob(iglob,jglob) + elem%matm (e,n1,n2)
              //m_glob.Set(iglob,jglob,m_glob.getVal(iglob,jglob)+ m_ematm[e*m_nodxelem + n1*m_nodxelem+n2]);
              
              m_glob.Set(iglob,jglob,(m_glob.getVal(iglob,jglob)+ tempm.getVal(n1,n2))*f*rho);
        }
      }
    } //FOR e //NOT PARALLEL
    
    double *mdiag_h = new double[m_node_count];

    //printf("mdiag: %f\n");
    for(int n =0; n<m_nodxelem;n++) {
      mdiag_h[n] = 0.0;
      for(int ig =0; ig<m_nodxelem;ig++){       
        mdiag_h[n] += m_glob.getVal(ig,n);
      }
      //printf("%f ",mdiag_h[n]); 
    }
    //printf("\n");
    
    memcpy_t(this->m_mdiag, mdiag_h, sizeof(double) * m_node_count);   
    
    delete mdiag_h;
  }
  double totmass = 0.0;
  /// DEVICE FUNCTION
  //// THIS ASSEMBLES LOOPING THROUGH NODES INSTEAD OF ELEMENTS
  //  CALCULATE ALSO M_DIAG
  dev_t void Domain_d::assemblyMassMatrix(){
  par_loop(n, m_node_count){
      
      double diag = 0.0;      
      for (int e=0; e<m_nodel_count[n];e++) {
        int eglob   = m_nodel     [m_nodel_offset[n]+e];
        int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        int offset  = eglob * m_nodxelem * m_nodxelem;
         //printf("glob n1, loc ne\n",n,ne);
        
        for (int n2=0;n2<m_node_count;n2++){
          for (int e2=0; e2<m_nodel_count[n2];e2++) {

            int eglob2   = m_nodel     [m_nodel_offset[n2]+e2];
            int ne2      = m_nodel_loc [m_nodel_offset[n2]+e2]; //LOCAL ELEMENT NODE INDEX m_nodel_local
            int offset2  = eglob2 * m_nodxelem * m_nodxelem;
            ////printf("glob ij: %d %d, loc ij %d, %d\n",n,n2,ne,ne2);
            double f = m_ematm[offset + ne * m_nodxelem + ne2];
            m_mglob[n*m_nodxelem + n2] += f;
            diag += f;
          }//nodel  count
          
        
        } // N2
        
        
      } // element
      m_mdiag [n] = diag;
      totmass += m_mdiag[n];
      //printf("m diag %d: %f\n",n, m_mdiag [n]);
    }//NODE N (N1)
  printf("TOTAL MASS: %f\n",totmass);
 }
  ///// INITIALIZE Displacements, velocities and Acceleration
  dev_t void Domain_d::InitUVA(){
    par_loop(n, m_node_count){    
      for (int d=0;d<m_dim;d++) {
        u[n*m_dim+d] = v[n*m_dim+d] = a[n*m_dim+d] = 0.0;
      }
    }
  }
  
//////////  CALCULATE MATM ///////////
///////// m_ematm //////
////  SEE BENSON 2017
//// Explicit Finite Element Methods for Large
//// Deformation Problems in Solid Mechanics
//// Section 
dev_t void Domain_d::calcElemMassMat() {   
  // !!!!! PREVIOUSLY JACOBIAN DETERMINANT SHOULD BE CALCULATED
// subroutine calculate_element_shapeMat ()
  // integer :: e,d
  // ! !rg=gauss[ig]
  // ! !sg=gauss[jg]
  // real(fp_kind), dimension(dim,nodxelem) :: dHrs
  // real(fp_kind), dimension(nodxelem,dim) :: x2
  // real(fp_kind), dimension(dim,dim) :: test
  // real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  // integer :: i,j,k, gp
  // real(fp_kind):: r, w, f
  // real(fp_kind), dimension(8,3):: gpc !!! gauss point coordinates, r,s,t
  // !! Update x2 vector (this is useful for strain and stress things)
  
  // r = 1.0/sqrt(3.0)
  // elem%matm(e,:,:) = 0.0d0
  // COPY THE DETERMINANT FROM DEV TO HOST
  par_loop (e, m_elem_count) {

  //Matrix *temph = new Matrix(m_dim, m_nodxelem);
  Matrix temph(m_dim, m_nodxelem);
  Matrix tempm(m_nodxelem, m_nodxelem);
  //Matrix *tempm = new Matrix(m_nodxelem, m_nodxelem);

  int offset = e * m_nodxelem * m_nodxelem;
  

  if (m_gp_count == 1 ) {  
  // do e=1, elem_count
    // gp = 1
    // if (elem%gausspc(e) .eq. 1) then
      // w = 2.0d0*dim
      // f = 1.0d0/(2.0d0*dim)
      // elem%math(e,gp,:,:) = 0.0d0
      // !!if (dim .eq. 2) then 
      
      // !!else !!!DIM 3
          // temph(:,:) = 0.0d0
          // do k =1, nodxelem
              // !temph(d,dim*(k-1)+d) = 0.125 !TODO: CHANGE IN 3D
              // elem%math(e,gp,:,k)=f
              // !print *, "temp h i ", d, "j ", dim*(k-1)+d, ": "
          // end do
          // ! print *, "temp h", temph(:,:)
      double w = pow(2.0, m_dim);
      double f = 1.0/pow(2.0, m_dim)*m_detJ[e] * rho[e];
      for (int k=0;k<m_nodxelem;k++){  
        for (int d=0;d<m_dim;d++)
          temph.Set(d,k,f);
          //printf("Element %d node %d mass: %f\n",e,k,f);
      }
      // TODO: DO NOT PERFORM THIS MULT
      
      //MatMul(temph.getTranspose(),temph,&tempm);
      //MatMul(temph->Transpose(),*temph,tempm);
      for (int i=0;i<m_nodxelem;i++)
        for (int j=0;j<m_nodxelem;j++)
          tempm.Set(i,j,f);
      ////printf("print hmat\n");
      //temph.Print();
      ////printf("print transpose\n");
      //temph.getTranspose().Print();
      ////printf("mmat \n");
      //temph->Transpose().Print();
      //tempm.Print();
					// //*jacob = 0.125 * MatMul(*dHrs,*x2);
          // MatMul(*dHrs,*x2,jacob);
          
          // ! temph(1,:) = [1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0]
          // ! temph(2,:) = [0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0]
          // ! temph(3,:) = [0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0]

      // !!end if  !!!!dim
      for (int i=0;i<m_nodxelem;i++)
        for (int j=0;j<m_nodxelem;j++){
          ////printf("rho %f\n", rho[e]);
          m_ematm[offset + i * m_nodxelem + j]=tempm.getVal(i,j); /*tempm.getVal(i,j) * m_detJ[e] * rho[e]*/
          ////printf("detJ %e ", tempm.getVal(i,j) * m_detJ[e]);
        }
        // elem%matm(e,:,:) = matmul(transpose(elem%math(e,gp,:,:)),elem%math(e,gp,:,:))*elem%rho(e,gp)*elem%detJ(e,gp)*w !!!2.0 ^3 WEIGHT
        
        // !print *, "MAT M", elem%matm(e,:,:)
      
    // else !!!!! GP > 1
      // if (dim .eq. 2) then 
        // w = 1.
        // else !!!DIM 3
          // gpc(1,:)=[-r,-r,-r];   gpc(2,:)=[ r,-r,-r];      gpc(3,:)=[-r, r,-r];      gpc(4,:)=[ r, r,-r];
          // gpc(5,:)=[-r,-r, r];   gpc(6,:)=[ r,-r, r];      gpc(7,:)=[-r, r, r];      gpc(8,:)=[ r, r, r];
          // do gp = 1,8
            // elem%math(e,gp, 1,:) = 0.125*[(1-gpc(gp,1))*(1-gpc(gp,2))*(1-gpc(gp,3)),(1+gpc(gp,1))*(1-gpc(gp,2))*(1-gpc(gp,3)), &
                                // (1+gpc(gp,1))*(1+gpc(gp,2))*(1-gpc(gp,3)),(1-gpc(gp,1))*(1+gpc(gp,2))*(1-gpc(gp,3)), &
                                // (1-gpc(gp,1))*(1-gpc(gp,2))*(1+gpc(gp,3)),(1+gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3)), &
                                // (1+gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3)),(1-gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3))]
            // !print *, "gp ",gp,  "detJ(e,gp)", elem%detJ(e,gp)
            // elem%matm(e,:,:) = elem%matm(e,:,:) + &
                               // matmul(transpose(elem%math(e,gp,:,:)),elem%math(e,gp,:,:))*elem%rho(e,gp)*elem%detJ(e,gp)*w !!!2.0 ^3 WEIGHT
          // end do
            // ! k = 1
            // ! do while (k <= nodxelem)
              // ! temph(2,2*k) = temph(1,2*k-1) !TODO: CHANGE IN 3D
              // ! k = k + 1
            // ! end do
 
          // ! elem%math(e,gp,:,:) = temph(:,:)*elem%detJ(e,gp)
          // !print *, "mat h ", elem%math(e,gp,:,:)        
      // end if !!! if dim
    // end if ! GP>1
  // end do !element
// end subroutine

  }//if gp == 1
  
    //delete temph;
    //delete tempm;
  }//par loop element

}

__global__ void calcElemMassMatKernel(Domain_d *dom_d){
  dom_d->calcElemMassMat();
}

__global__ void assemblyMassMatrixKernel(Domain_d *dom_d){
  dom_d->assemblyMassMatrix();
}

}//MetFEM
