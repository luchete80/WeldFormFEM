#include "NodeData.cuh"


 void AllocateNodeData(NodeData *node, const int &nc){
  //node_count = nc;
  cudaMalloc((void **)&node->x, nc * sizeof (double3));
  cudaMalloc((void **)&node->v, nc * sizeof (double3));
  cudaMalloc((void **)&node->a, nc * sizeof (double3));
  cudaMalloc((void **)&node->u, nc * sizeof (double3));
 }