#include "ElementData.cuh"


//// FOR A SINGLE ELEMENT DATA

void AllocateElementData(ElementData *elem, const int &dim, const int &el_count, const int &gp, const int &nodxelem){
  
  cudaMalloc((void **)&elem->pressure, el_count * sizeof (double)); //8 values per dim 
  cudaMalloc((void **)&elem->gausspc, el_count * sizeof (int)); 
  
  cudaMalloc((void **)&elem->dHxy_detJ, el_count*gp*dim*nodxelem*sizeof (double));

  //ATTENTION; THIS SHOULD BE CHANGE WHEN ADDING DIFFERENT TYPEOFELEMENTS
  cudaMalloc((void **)&elem->elnod, el_count * nodxelem * sizeof (unsigned long));
  cudaMalloc((void **)&elem->elnod_offset, el_count * sizeof (unsigned long));
  
  
  //Assuming reduced integration
  //cudaMalloc((void **)&elem->dHxy_detJ, el_count * dim * nodxelem* sizeof (double)); //8 values per dim 


 
     // cudaMalloc((void **)&elem->elnod(el_count,nodxelem))

    // cudaMalloc((void **)&elem->dof(el_count,dim*nodxelem))
    // cudaMalloc((void **)&elem->vol(el_count))
    // cudaMalloc((void **)&elem->vol_inc(el_count))
    // cudaMalloc((void **)&elem->vol_0(el_count))
    // cudaMalloc((void **)&elem->x2(el_count,nodxelem,dim))
    // cudaMalloc((void **)&elem->jacob(el_count,gp,dim,dim))
    // cudaMalloc((void **)&elem->detj(el_count,gp))
    // cudaMalloc((void **)&elem->sigma_eq(el_count,gp)) !But is constant??
    // cudaMalloc((void **)&elem->dHxy(el_count,gp,dim,nodxelem))
    // cudaMalloc((void **)&elem->dHxy_detJ(el_count,gp,dim,nodxelem)) !!!! STORE LIKE THIS TO SAVE CALCULATION TIME (THIS IS USED  TO CALC FORCES INTEGRATING IT )
    // cudaMalloc((void **)&elem->dHxy0(el_count,gp,dim,nodxelem)) !!!USED FOR DEFORMATION GRADIENT ONLY FOR FULL INTEGRATION ELEMENTS 
    // cudaMalloc((void **)&elem->dHrs(el_count,gp,dim,nodxelem))
    // cudaMalloc((void **)&elem->sigma(el_count,gp,dim,dim))  !!!THIS IS A DIMxDIM SYMMETRIC TENSOR

    // cudaMalloc((void **)&elem->uele (el_count,dim*nodxelem,1)) 

    // cudaMalloc((void **)&elem->vele (el_count,dim*nodxelem,1)) 
    
    // cudaMalloc((void **)&elem->mass(el_count)) !Mass matrix    
    
    // cudaMalloc((void **)&elem->c_s(el_count,gp))
    // cudaMalloc((void **)&elem->p_visc(el_count,gp))
    // cudaMalloc((void **)&elem->e_length(el_count))

    // cudaMalloc((void **)&elem->matm(el_count,nodxelem,nodxelem)) !Mass matrix
    // cudaMalloc((void **)&elem->math(el_count,gp,1,nodxelem)) !Mass matrix
    
    // cudaMalloc((void **)&elem->hourg_nodf(el_count,nodxelem,dim)) !AS 1 COLUMN OR NOT????? Mass matrix
    
    // cudaMalloc((void **)&elem->f_int(el_count,nodxelem,dim))
    // cudaMalloc((void **)&elem->f_ext(el_count,nodxelem,dim))
    
    // cudaMalloc((void **)&elem->rho(el_count,gp)) !AT FIRST ONLY ONE POINT
    // cudaMalloc((void **)&elem->rho_0(el_count,gp))
    // cudaMalloc((void **)&elem->pressure(el_count,gp))
    // cudaMalloc((void **)&elem->cs(el_count))
    // cudaMalloc((void **)&elem->shear_stress(el_count,gp, dim,dim))
    // cudaMalloc((void **)&elem->str_rate(el_count,gp, dim,dim))
    // cudaMalloc((void **)&elem->str_inc(el_count,gp, dim,dim))
    // cudaMalloc((void **)&elem->rot_rate(el_count,gp, dim,dim))
      
    // if (Dim .eq. 2) then
      // cudaMalloc((void **)&elem->bl (el_count,gp,3,dim*nodxelem))
      // cudaMalloc((void **)&elem->bnl(el_count,gp, 4,dim*nodxelem))
      // cudaMalloc((void **)&elem->strain(el_count,gp, 4,1))
      // !cudaMalloc((void **)&elem->str_rate(el_count,gp, 4,1))
      // !cudaMalloc((void **)&elem->rot_rate(el_count,gp, 4,1))
    // else 
      // cudaMalloc((void **)&elem->bl (el_count,gp,6,dim*nodxelem)) 
      // cudaMalloc((void **)&elem->strain(el_count,gp, 6,1)) !!VECTORIZED 
      // !cudaMalloc((void **)&elem->str_rate(el_count,gp, 6,1))
      // !cudaMalloc((void **)&elem->rot_rate(el_count,gp, 6,1))
    // end if 
    
    // elem->gausspc(:) = gp
    
}