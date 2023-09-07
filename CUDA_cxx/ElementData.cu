#include "ElementData.cuh"


void AllocateElements(ElementData *elem, const int &dim, const int &el_count){
  
  cudaMalloc((void **)&elem->pressure, el_count * sizeof (double)); //8 values per dim 
  
  //Assuming reduced integration
  //cudaMalloc((void **)&elem->dHxy_detJ, el_count * dim * nodxelem* sizeof (double)); //8 values per dim 
}