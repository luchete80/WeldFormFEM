#include "Domain_d.h"

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
          m_f[n*m_dim + d] += m_f_elem[offset + ne*m_dim + d];
        }
      }
    } // element

  }//assemblyForcesNonLock


}//MetFEM