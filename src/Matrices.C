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


////HERE THE LOOP IS PER NODE INSTEAD THAN PER ELEMENT
  dev_t void Domain_d::assemblyForces(){

    //if ()
  par_loop(n, m_node_count){
      for (int e=0; e<m_nodel_count[n];e++) {
        int eglob   = m_nodel     [m_nodel_offset[n]+e];
        int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        int offset  = e * m_nodxelem * m_dim;
        for (int d=0;d<m_dim;d++){
          //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
          m_fi[n*m_dim + d] += m_f_elem[offset + ne*m_dim + d];
        }
      }
    } // element

  }//assemblyForcesNonLock


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

  host_ Domain_d::calcMassDiagFromElementNodes(int *elnod) // To use existing array
  {
    Matrix m_glob(m_node_count,m_node_count);
    int n1,n2, iglob,jglob;
    // m_glob (:,:) = 0.0d0
      for (int e=0;e<m_elem_count;e++){
        // !print *, "elem ", e
        for( n1 =0; n1<m_nodxelem;n1++){
          for( n2 =0; n2<m_nodxelem;n2++){
                // !print *, "elem ", e, "node ", n, " i j matm ",i, j, elem%matm (e,dim*(n-1)+i,dim*(n2-1)+j)            
                // iglob  = elnod(e,n1) 
                // jglob  = elnod(e,n2) 
                // ! print *, "iloc, jloc ",dim*(n-1)+i, dim*(n2-1)+j, "iglob, jglob", iglob,jglob
                // m_glob(iglob,jglob) = m_glob(iglob,jglob) + elem%matm (e,n1,n2)
                //m_glob->Set(iglob,jglob,m_glob->getVal(iglob,jglob)+;
          }
        }
      }
  }
  
  ///// INITIALIZE Displacements, velocities and Acceleration
  dev_t void Domain_d::InitUVA(){
    par_loop(n, m_node_count){    
      for (int d=0;d<m_dim;d++) {
        u[n*m_dim+d] = v[n*m_dim+d] = a[n*m_dim+d] = 0.0;
      }
    }
  }

}//MetFEM